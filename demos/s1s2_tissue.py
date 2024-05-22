# # S1S2 protocol on tissue
#

from typing import NamedTuple
from pathlib import Path
import beat.single_cell
import dolfin
import matplotlib.pyplot as plt
import numpy as np
import pint
import gotranx
import beat

try:
    import ufl_legacy as ufl
except ImportError:
    import ufl


ureg = pint.UnitRegistry()


class Geometry(NamedTuple):
    mesh: dolfin.Mesh
    f0: dolfin.Constant
    s0: dolfin.Constant


def define_conductivity_tensor(
    chi,
    f0,
    s0,
    g_il=0.17,
    g_it=0.019,
    g_el=0.62,
    g_et=0.24,
):
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = g_il * ureg("S/m")
    sigma_it = g_it * ureg("S/m")
    sigma_el = g_el * ureg("S/m")
    sigma_et = g_et * ureg("S/m")

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(chi)
    s_l = (sigma_l / chi).to("uA/mV").magnitude
    s_t = (sigma_t / chi).to("uA/mV").magnitude

    # Define conductivity tensor

    A = dolfin.as_matrix(
        [
            [f0[0], s0[0]],
            [f0[1], s0[1]],
        ],
    )

    M_star = ufl.diag(dolfin.as_vector([s_l, s_t]))

    M = A * M_star * A.T

    return M


def setup_geometry(dx, Lx, Ly):

    mesh = dolfin.RectangleMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0),
        dolfin.Point(Lx, Ly),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
    )

    return mesh


def load_timesteps_from_xdmf(xdmffile):
    import xml.etree.ElementTree as ET

    times = {}
    i = 0
    tree = ET.parse(xdmffile)
    for elem in tree.iter():
        if elem.tag == "Time":
            times[i] = float(elem.get("Value"))
            i += 1

    return times


def delete_old_file(file: str):
    Path(file).unlink(missing_ok=True)
    Path(file).with_suffix(".h5").unlink(missing_ok=True)


def run_model(
    data: Geometry,
    resultsdir: Path,
    num_beats: int = 5,
    dt: float = 0.05,
    save_freq: int = 20,
    g_il=0.16069,
    g_el=0.625,
    g_it=0.04258,
    g_et=0.236,
    **kwargs,
):
    print("Running model")
    # Load the model
    model_path = Path("tentusscher_panfilov_2006_epi_cell.py")
    if not model_path.is_file():
        here = Path(__file__).parent
        ode = gotranx.load_ode(
            here
            / ".."
            / "odes"
            / "tentusscher_panfilov_2006"
            / "tentusscher_panfilov_2006_epi_cell.ode"
        )
        code = gotranx.cli.gotran2py.get_code(
            ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
        )
        model_path.write_text(code)

    import tentusscher_panfilov_2006_epi_cell

    model = tentusscher_panfilov_2006_epi_cell.__dict__

    # Surface to volume ratio
    chi = 1400.0 * ureg("cm**-1")

    # Membrane capacitance
    C_m = 1.0 * ureg("uF/cm**2")

    fun = model["forward_generalized_rush_larsen"]

    y = model["init_state_values"]()

    time = dolfin.Constant(0.0)
    CI1 = 350.0
    CI0 = 50.0
    CI = dolfin.Constant(CI1)
    CIincr = 50.0
    parameters = model["init_parameter_values"](stim_amplitude=0.0)

    I_s_expr = dolfin.Expression(
        "((std::fmod(time,PCL) >= start) & (std::fmod(time,PCL) <= duration + start)) ? amplitude : 0.0",
        time=time,
        start=5.0,
        duration=5.0,
        amplitude=1.0,
        PCL=CI,
        degree=0,
    )
    subdomain_data = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_data.set_all(0)
    marker = 1
    dolfin.CompiledSubDomain("x[0] < 2 * dx", dx=dx).mark(subdomain_data, 1)
    ds = dolfin.Measure("ds", domain=mesh, subdomain_data=subdomain_data)(marker)
    I_s = beat.base_model.Stimulus(dz=ds, expr=I_s_expr)

    V_ode = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    parameters_ode = np.zeros((len(parameters), V_ode.dim()))
    parameters_ode.T[:] = parameters

    M = define_conductivity_tensor(
        chi=chi,
        f0=data.f0,
        s0=data.s0,
        g_il=g_il,
        g_it=g_it,
        g_el=g_el,
        g_et=g_et,
    )
    pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, C_m=C_m.magnitude)
    ode = beat.odesolver.DolfinODESolver(
        v_ode=dolfin.Function(V_ode),
        v_pde=pde.state,
        fun=fun,
        init_states=y,
        parameters=parameters_ode,
        num_states=len(y),
        v_index=model["state_index"]("V"),
    )
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

    fname = (resultsdir / "V.xdmf").as_posix()
    delete_old_file(fname)

    # Pick 12 uniformly distributed points on the mesh
    points = np.array(
        [
            [0.1, 0.1],
            [0.1, 0.5],
            [0.1, 0.9],
            [0.5, 0.1],
            [0.5, 0.5],
            [0.5, 0.9],
            [0.9, 0.1],
            [0.9, 0.5],
            [0.9, 0.9],
        ]
    )
    V_values: list[list[float]] = [[] for _ in points]
    times = []
    stim_values = []

    def save(t):
        v = solver.pde.state.vector().get_local()
        print(f"Solve for {t=:.2f}, {v.max() =}, {v.min() = }")
        with dolfin.XDMFFile(mesh.mpi_comm(), fname) as xdmf:
            xdmf.parameters["functions_share_mesh"] = True
            xdmf.parameters["rewrite_function_mesh"] = False
            xdmf.write_checkpoint(
                solver.pde.state,
                "V",
                float(t),
                dolfin.XDMFFile.Encoding.HDF5,
                True,
            )

        times.append(t)
        stim_values.append(I_s_expr(0))
        for i, p in enumerate(points):
            V_values[i].append(solver.pde.state(p))
        fig, ax = plt.subplots(2, 1)
        ax[0].plot(times, stim_values)
        for p in points:
            ax[1].plot(times, V_values[i], label=f"{p}")
        ax[1].set_xlabel("Time [ms]")
        ax[1].set_ylabel("V [mV]")
        ax[1].legend(title="Point")
        ax[0].set_title("Stimulus")
        fig.savefig(resultsdir / "V.png")
        plt.close()
        print("Saved to ", resultsdir / "V.png")
        np.save(
            resultsdir / "V.npy",
            {"V": V_values, "time": times, "stim": stim_values},
            allow_pickle=True,
        )

    t = 0.0
    save_freq = int(1.0 / dt)
    i = 0

    CIs = np.arange(CI1, CI0 - CIincr, -CIincr)

    for CI_value in CIs:
        CI.assign(CI_value)
        for _ in range(num_beats):
            for ti in np.arange(0, CI_value, dt):
                if i % save_freq == 0:
                    save(t)

                solver.step((t, t + dt))
                i += 1
                t += dt

                time.assign(ti)


def get_microstructure(
    dim: int, transverse: bool = False
) -> tuple[dolfin.Constant, dolfin.Constant]:

    if transverse:
        f0 = dolfin.Constant((0.0, 1.0))
        s0 = dolfin.Constant((1.0, 0.0))
    else:
        f0 = dolfin.Constant((1.0, 0.0))
        s0 = dolfin.Constant((0.0, 1.0))

    return f0, s0


results_folder = Path("results-s1s2-tissue")
save_every_ms = 1.0
dimension = 2
transverse = False
end_time = 20.0
dt = 0.05
overwrite = False
stim_amp = 5000.0
mesh_unit = "cm"

# Load mesh
mesh_unit = mesh_unit

dx = 0.02 * ureg("cm").to(mesh_unit).magnitude
L = 1.0 * ureg("cm").to(mesh_unit).magnitude
mesh = setup_geometry(Lx=L, Ly=L, dx=dx)


V = dolfin.FunctionSpace(mesh, "CG", 1)


g_il = 0.16069
g_el = 0.625
g_it = 0.04258
g_et = 0.236

f0, s0 = get_microstructure(dimension, transverse)


data = Geometry(
    mesh=mesh,
    f0=f0,
    s0=s0,
)
save_freq = round(save_every_ms / dt)
run_model(
    data=data,
    g_il=g_il,
    g_it=g_it,
    g_el=g_el,
    g_et=g_et,
    save_freq=save_freq,
    resultsdir=results_folder,
    dt=dt,
    stim_amp=stim_amp,
    mesh_unit=mesh_unit,
)
