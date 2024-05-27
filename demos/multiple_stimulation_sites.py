# # Multiple stimulation sites
#

from pathlib import Path

import dolfin
import matplotlib.pyplot as plt
import numpy as np

import gotranx
import beat
import pyvista
import beat.viz


class StimulationProtocol(dolfin.UserExpression):
    def __init__(self, time, duration=5.0, delay=50.0, **kwargs):
        self.time = time
        self.duration = duration
        self.delay = delay
        super().__init__(**kwargs)

    def eval_cell(self, value, x, ufc_cell):
        t = float(self.time)
        if t < self.duration and x[0] <= 0.05:
            value[0] = 20.0

        elif (
            t > self.delay
            and t < self.delay + self.duration
            and (x[0] - 0.8) ** 2 + (x[1] - 0.8) ** 2 <= 0.05**2
        ):
            value[0] = 40.0
        else:
            value[0] = 0.0


def setup_geometry(dx, Lx, Ly):

    mesh = dolfin.RectangleMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0),
        dolfin.Point(Lx, Ly),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
    )

    return mesh


def run_model(
    data: beat.Geometry,
    resultsdir: Path,
    num_beats: int = 5,
    dt: float = 0.05,
    save_freq: int = 20,
    g_il=0.16069,
    g_el=0.625,
    g_it=0.04258,
    g_et=0.236,
    end_time: float = 5000.0,
    **kwargs,
):
    print("Running model")
    # Load the model
    model_path = Path("courtemanche_ramirez_nattel_1998.py")
    if not model_path.is_file():
        here = Path(__file__).parent
        ode = gotranx.load_ode(
            here / ".." / "odes" / "courtemanche_ramirez_nattel_1998.ode"
        )
        code = gotranx.cli.gotran2py.get_code(
            ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
        )
        model_path.write_text(code)

    import courtemanche_ramirez_nattel_1998

    model = courtemanche_ramirez_nattel_1998.__dict__

    # Surface to volume ratio
    chi = 1400.0 * beat.units.ureg("cm**-1")

    # Membrane capacitance
    C_m = 1.0 * beat.units.ureg("uF/cm**2")

    fun = model["forward_generalized_rush_larsen"]

    y = model["init_state_values"]()

    time = dolfin.Constant(0.0)
    parameters = model["init_parameter_values"](stim_amplitude=0.0)

    I_s = beat.base_model.Stimulus(dz=dolfin.dx, expr=StimulationProtocol(time))

    V_ode = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    parameters_ode = np.zeros((len(parameters), V_ode.dim()))
    parameters_ode.T[:] = parameters

    M = beat.conductivities.define_conductivity_tensor(
        chi=chi,
        f0=data.f0,
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
    beat.postprocess.delete_old_file(fname)

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
        for i, p in enumerate(points):
            V_values[i].append(solver.pde.state(p))
        fig, ax = plt.subplots()
        for p in points:
            ax.plot(times, V_values[i], label=f"{p}")
        ax.set_xlabel("Time [ms]")
        ax.set_ylabel("V [mV]")
        ax.legend(title="Point")
        fig.savefig(resultsdir / "V.png")
        plt.close()
        print("Saved to ", resultsdir / "V.png")
        np.save(
            resultsdir / "V.npy",
            {"V": V_values, "time": times},
            allow_pickle=True,
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
    transverse: bool = False,
) -> tuple[dolfin.Constant, dolfin.Constant]:

    if transverse:
        f0 = dolfin.Constant((0.0, 1.0))
        s0 = dolfin.Constant((1.0, 0.0))
    else:
        f0 = dolfin.Constant((1.0, 0.0))
        s0 = dolfin.Constant((0.0, 1.0))

    return f0, s0


results_folder = Path("results-multiple-stimulation-sites")
save_every_ms = 1.0
transverse = False
end_time = 20.0
dt = 0.05
overwrite = False
stim_amp = 5000.0
mesh_unit = "cm"

# Load mesh
mesh_unit = mesh_unit

dx = 0.05 * beat.units.ureg("cm").to(mesh_unit).magnitude
L = 1.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
mesh = setup_geometry(Lx=L, Ly=L, dx=dx)


V = dolfin.FunctionSpace(mesh, "CG", 1)


g_il = 0.16069
g_el = 0.625
g_it = 0.04258
g_et = 0.236

f0, s0 = get_microstructure(transverse)


data = beat.Geometry(
    mesh=mesh,
    f0=f0,
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
    end_time=end_time,
)