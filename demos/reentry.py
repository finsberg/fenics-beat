# # Conduction velocity and ECG for slabs
# In this demo we will show how to compute conduction velocity and ECG for a Slab geometry.
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
import pyvista
import beat.viz
from beat.single_cell import get_steady_state

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
    V_dg,
    indices,
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

    s_l_func = dolfin.Function(V_dg)
    s_l_func.vector()[:] = s_l
    s_l_func.vector()[indices] = 1e-7
    s_t_func = dolfin.Function(V_dg)
    s_t_func.vector()[:] = s_t
    s_t_func.vector()[indices] = 1e-7

    # Define conductivity tensor

    A = dolfin.as_matrix(
        [
            [f0[0], s0[0]],
            [f0[1], s0[1]],
        ],
    )

    M_star = ufl.diag(dolfin.as_vector([s_l_func, s_t_func]))

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
    end_time: float = 15.0,
    dt: float = 0.05,
    save_freq: int = 20,
    g_il=0.16069,
    g_el=0.625,
    g_it=0.04258,
    g_et=0.236,
    **kwargs,
):
    print("Running model")
    # # Load the model
    model_path = Path("courtemanche_ramirez_nattel_1998.py")
    # model_path = Path("tentusscher_panfilov_2006_epi_cell.py")
    if not model_path.is_file():
        here = Path(__file__).parent
        ode = gotranx.load_ode(
            here / ".." / "odes" / "courtemanche_ramirez_nattel_1998.ode"
        )
        #     here = Path(__file__).parent
        #     ode = gotranx.load_ode(
        #         here
        #         / ".."
        #         / "odes"
        #         / "tentusscher_panfilov_2006"
        #         / "tentusscher_panfilov_2006_epi_cell.ode"
        #     )
        code = gotranx.cli.gotran2py.get_code(
            ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
        )
        model_path.write_text(code)

    # import tentusscher_panfilov_2006_epi_cell

    # model = tentusscher_panfilov_2006_epi_cell.__dict__

    import courtemanche_ramirez_nattel_1998

    model = courtemanche_ramirez_nattel_1998.__dict__

    # Surface to volume ratio
    chi = 1400.0 * ureg("cm**-1")

    # Membrane capacitance
    C_m = 1.0 * ureg("uF/cm**2")

    fun = model["forward_generalized_rush_larsen"]

    y = model["init_state_values"]()

    time = dolfin.Constant(0.0)
    parameters = model["init_parameter_values"](stim_amplitude=0.0)
    # parameters = model["init_parameter_values"](amp=0.0)

    I_s_expr = dolfin.Expression(
        "((std::fmod(time,PCL) >= start) & (std::fmod(time,PCL) <= duration + start)) ? amplitude : 0.0",
        time=time,
        start=5.0,
        duration=5.0,
        amplitude=1.0,
        PCL=150.0,
        degree=0,
    )
    subdomain_data = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_data.set_all(0)
    marker = 1
    dolfin.CompiledSubDomain(
        "std::pow((x[0] - 0.2), 2) + std::pow((x[1] - 0.5), 2) < 0.05", dx=dx
    ).mark(subdomain_data, 1)

    ds = dolfin.Measure("ds", domain=mesh, subdomain_data=subdomain_data)(marker)
    I_s = beat.base_model.Stimulus(dz=ds, expr=I_s_expr)

    V_ode = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    parameters_ode = np.zeros((len(parameters), V_ode.dim()))
    parameters_ode.T[:] = parameters

    fibrosis = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    fibrosis.set_all(0)
    marker = 1
    dolfin.CompiledSubDomain(
        "std::pow((x[0] - 0.5),2) + std::pow((x[1] - 0.5), 2) < 0.05", dx=dx
    ).mark(fibrosis, marker)

    indices = np.where(fibrosis.array() == 1)[0]
    np.random.shuffle(indices)
    N = len(indices)
    indices1 = indices[: int(N * 0.7)]
    indices2 = indices[int(N * 0.7) :]

    V_dg = dolfin.FunctionSpace(mesh, "DG", 0)

    g_K1_index = model["parameter_index"]("g_K1")
    # g_K1_index = model["parameter_index"]("scale_drug_IK1")
    g_K1_value = parameters[g_K1_index]
    g_K1_func = dolfin.Function(V_dg)
    g_K1_func.vector()[:] = g_K1_value
    g_K1_func.vector()[indices1] = g_K1_value * 0.5
    parameters_ode[g_K1_index, :] = (
        dolfin.interpolate(g_K1_func, V_ode).vector().get_local()
    )

    g_Na_index = model["parameter_index"]("g_Na")
    # g_Na_index = model["parameter_index"]("scale_drug_INa")
    g_Na_value = parameters[g_Na_index]
    g_Na_func = dolfin.Function(V_dg)
    g_Na_func.vector()[:] = g_Na_value
    g_Na_func.vector()[indices1] = g_Na_value * 0.6
    parameters_ode[g_Na_index, :] = (
        dolfin.interpolate(g_Na_func, V_ode).vector().get_local()
    )

    g_CaL_index = model["parameter_index"]("g_Ca_L")
    # g_CaL_index = model["parameter_index"]("scale_drug_ICaL")
    g_CaL_value = parameters[g_CaL_index]
    g_CaL_func = dolfin.Function(V_dg)
    g_CaL_func.vector()[:] = g_CaL_value
    g_CaL_func.vector()[indices1] = g_CaL_value * 0.5
    parameters_ode[g_CaL_index, :] = (
        dolfin.interpolate(g_CaL_func, V_ode).vector().get_local()
    )

    M = define_conductivity_tensor(
        chi=chi,
        V_dg=V_dg,
        indices=indices2,
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

    t = 0.0
    save_freq = int(1.0 / dt)

    i = 0
    while t < end_time + 1e-12:
        # Make sure to save at the same time steps that is used by Ambit

        if i % save_freq == 0:
            save(t)

        solver.step((t, t + dt))
        i += 1
        t += dt


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


results_folder = Path("results-reentry")
save_every_ms = 1.0
dimension = 2
transverse = False
end_time = 200.0
dt = 0.05
overwrite = False
stim_amp = 5000.0
mesh_unit = "cm"

# Load mesh
mesh_unit = mesh_unit

dx = 0.01 * ureg("cm").to(mesh_unit).magnitude
L = 1.0 * ureg("cm").to(mesh_unit).magnitude
mesh = setup_geometry(Lx=L, Ly=L, dx=dx)


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
    end_time=end_time,
    dt=dt,
    stim_amp=stim_amp,
    mesh_unit=mesh_unit,
)
