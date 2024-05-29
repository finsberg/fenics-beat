# # S1S2 protocol on tissue
#

from typing import NamedTuple
from pathlib import Path
import beat.single_cell
import dolfin
import matplotlib.pyplot as plt
import numpy as np
import gotranx
import beat


results_folder = Path("results-s1s2-tissue")
save_every_ms = 1.0
dimension = 2
transverse = False
end_time = 20.0
dt = 0.05
overwrite = False
stim_amp = 500.0
mesh_unit = "cm"

# Load mesh
mesh_unit = mesh_unit

dx = 0.02 * beat.units.ureg("cm").to(mesh_unit).magnitude
L = 1.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
data = beat.geometry.get_2D_slab_geometry(Lx=L, Ly=L, dx=dx)


V = dolfin.FunctionSpace(data.mesh, "CG", 1)

save_freq = round(save_every_ms / dt)

print("Running model")

# Load the model
model_path = Path("tentusscher_panfilov_2006_epi_cell.py")
if not model_path.is_file():
    here = Path.cwd()
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
chi = 1400.0 * beat.units.ureg("cm**-1")

# Membrane capacitance
C_m = 1.0 * beat.units.ureg("uF/cm**2")

fun = model["forward_generalized_rush_larsen"]

y = model["init_state_values"]()

time = dolfin.Constant(0.0)
CI1 = 350.0
CI0 = 50.0
CI = dolfin.Constant(CI1)
CIincr = 50.0
parameters = model["init_parameter_values"](stim_amplitude=0.0)


subdomain_data = dolfin.MeshFunction(
    "size_t", data.mesh, data.mesh.topology().dim() - 1
)
subdomain_data.set_all(0)
marker = 1
dolfin.CompiledSubDomain("x[0] < 2 * dx", dx=dx).mark(subdomain_data, 1)


I_s = beat.stimulation.define_stimulus(
    mesh=data.mesh,
    chi=chi,
    time=time,
    subdomain_data=subdomain_data,
    marker=marker,
    amplitude=stim_amp,
    mesh_unit=mesh_unit,
    PCL=CI,
    start=5.0,
    duration=5.0,
)

V_ode = dolfin.FunctionSpace(data.mesh, "Lagrange", 1)
parameters_ode = np.zeros((len(parameters), V_ode.dim()))
parameters_ode.T[:] = parameters

M = beat.conductivities.define_conductivity_tensor(
    f0=data.f0,
    **beat.conductivities.default_conductivities("Niederer"),
)
pde = beat.MonodomainModel(time=time, mesh=data.mesh, M=M, I_s=I_s, C_m=C_m.magnitude)
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

fname = (results_folder / "V.xdmf").as_posix()
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

    stim_values.append(I_s.expr(0))
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
    fig.savefig(results_folder / "V.png")
    plt.close()
    print("Saved to ", results_folder / "V.png")
    np.save(
        results_folder / "V.npy",
        {"V": V_values, "time": times, "stim": stim_values},
        allow_pickle=True,
    )


save_freq = int(1.0 / dt)
num_beats = 5

CIs = np.arange(CI1, CI0 - CIincr, -CIincr)


def run():
    t = 0.0
    i = 0
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
                if t > end_time:
                    print("Terminating")
                    return


run()
