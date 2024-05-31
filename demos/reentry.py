# # Inducing reentry in a 2D sheet of cardiac tissue
#
# In this demo we show how to induce reentry in a 2D sheet of cardiac tissue using two stimulations, i.e one wave where we stimulate the left side of the square and then one second stimulation where we stimulate the third quadrant of the square.
#
# First we make the necessary imports

from pathlib import Path
import dolfin
import matplotlib.pyplot as plt
import numpy as np
import gotranx
import beat

# And set up some basic parameters

results_folder = Path("results-reentry")
save_every_ms = 1.0
transverse = False
end_time = 3000.0
dt = 0.05
save_freq = round(save_every_ms / dt)
overwrite = False
stim_amp = 5000.0
mesh_unit = "cm"
dx = 0.4 * beat.units.ureg("mm").to(mesh_unit).magnitude
L = 5.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
data = beat.geometry.get_2D_slab_geometry(
    Lx=L,
    Ly=L,
    dx=dx,
    transverse=transverse,
)

V = dolfin.FunctionSpace(data.mesh, "CG", 1)

# For this simulation we will use a model from {courtemanche1998ionic}`courtemanche1998ionic`. The model is taken form https://models.physiomeproject.org/workspace/courtemanche_ramirez_nattel_1998
#

# +
print("Running model")
# Load the model
model_path = Path("courtemanche_ramirez_nattel_1998.py")
if not model_path.is_file():
    here = Path.cwd()
    ode = gotranx.load_ode(
        here / ".." / "odes" / "courtemanche_ramirez_nattel_1998.ode"
    )
    code = gotranx.cli.gotran2py.get_code(
        ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
    )
    model_path.write_text(code)

import courtemanche_ramirez_nattel_1998

model = courtemanche_ramirez_nattel_1998.__dict__
# -


# Membrane capacitance
C_m = 1.0 * beat.units.ureg("uF/cm**2")

# Here we also reduce the conductance for the sodium channel to make sure that the speed of the traveling wave is reduced
fun = model["forward_generalized_rush_larsen"]
y = model["init_state_values"]()
time = dolfin.Constant(0.0)
parameters = model["init_parameter_values"](stim_amplitude=0.0, g_Na=2.5)

# The duration of the stimulation is set to 5 ms and the second stimulation is done 400 ms after the initial stimulation. This timing is important to get spiral wave appearing
duration = 5.0
delay = 400.0
# We specify this stimulation protocol as an expressions

# +
expr = dolfin.Expression(
    "(t < duration && x[0] < 0.05) ? 20.0 : "
    "((t > delay && t < delay+duration) && (x[0] < L / 2 && x[1] < L / 2)) ? 20.0 : 0.0",
    t=time,
    L=L,
    delay=delay,
    duration=duration,
    degree=0,
)

I_s = beat.base_model.Stimulus(
    dz=dolfin.dx, expr=expr  # StimulationProtocol(time, L=L),
)
# -

# Next we set up the models

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

# and we also pick 12 uniformly distributed points on the mesh and keep track of the voltage in these points

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
    with dolfin.XDMFFile(data.mesh.mpi_comm(), fname) as xdmf:
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
    fig.savefig(results_folder / "V.png")
    plt.close()
    print("Saved to ", results_folder / "V.png")
    np.save(
        results_folder / "V.npy",
        {"V": V_values, "time": times},
        allow_pickle=True,
    )


# Finally we solve it

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


# ```{figure} ../docs/_static/reentry.mp4
# ---
# name: pvc
# ---
# ```

# ```{bibliography}
# ```
#
