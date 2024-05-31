# # Premature Ventricular Complexes (PVCs)
# In this demo we try to replicate the experimental setup conducted in {cite}`zhang2021mechanisms`
#
# First we do the necceary imports
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import beat
import dolfin
import gotranx
from beat.single_cell import get_steady_state

# Next we set the output directory for the results and define the geometry. Here we specify an interval mesh of 200 cells with 0.015 cm between each cell

here = Path.cwd()
outdir = here / "results-pvc"
mesh_unit = "cm"
dx = 0.015
num_cells = 200
L = num_cells * dx
mesh = dolfin.IntervalMesh(num_cells, 0, L)

# We will use the tensusscher panfilov model

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

# Here we can also specify whether the stimulation should be on one side or if we should stimulate all cell at once. We can also specify an end time for the simulation.

traveling_wave = False
# Change this to run the simulations for longer
end_time = 1500.0

# We specify a diffusion coefficient and a membrane capacitance

D = 0.0005 * beat.units.ureg("cm**2 / ms")
Cm = 1.0 * beat.units.ureg("uF/cm**2")

# Next we run a single cell model with the to get the correct steady state solutions.
# We run this for 50 beats with a stimulation every 500.0 ms
# +
parameters = model["init_parameter_values"](stim_period=500.0)

dt = 0.01
nbeats = 50
fun = model["forward_generalized_rush_larsen"]
y = get_steady_state(
    fun=fun,
    init_states=model["init_state_values"](),
    parameters=parameters,
    outdir=outdir / "prebeats",
    BCL=500,
    nbeats=nbeats,
    track_indices=[model["state_index"]("V"), model["state_index"]("Ca_i")],
    dt=dt,
)
# -

# Now we create the stimulus current which is 0.0 is we stimulate all the cells at once (in which case the stimulation is applied directly to all cells), or in the case of a traveling wave we stimulate the two first cells.
#
# In both cases we stimulate after 100 ms with period of 500 ms

time = dolfin.Constant(0.0)
if traveling_wave:
    parameters = model["init_parameter_values"](stim_amplitude=0.0)

    I_s_expr = dolfin.Expression(
        "std::fmod(time,PCL) >= start ? (std::fmod(time,PCL) <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=100.0,
        duration=2.0,
        amplitude=1.0,
        PCL=500.0,
        degree=0,
    )
    subdomain_data = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_data.set_all(0)
    marker = 1
    dolfin.CompiledSubDomain("x[0] < 2 * dx", dx=dx).mark(subdomain_data, 1)

    ds = dolfin.Measure("ds", domain=mesh, subdomain_data=subdomain_data)(marker)
    I_s = beat.base_model.Stimulus(dz=ds, expr=I_s_expr)
else:
    I_s = dolfin.Constant(0.0)
    parameters = model["init_parameter_values"](stim_start=100.0, stim_period=500.0)

# Now we would like to change the conductances for the cells in the right most part of the cable. We choose a first order Lagrange space for this, so that we have one set of parameters for each cells, and then we copy over the default parameters
V_ode = dolfin.FunctionSpace(mesh, "Lagrange", 1)
parameters_ode = np.zeros((len(parameters), V_ode.dim()))
parameters_ode.T[:] = parameters

# Next we set the values og $g_Kr$ to zero in the right part of the cable
g_Kr_index = model["parameter_index"]("g_Kr")
g_Kr_value = parameters[g_Kr_index]
parameters_ode[g_Kr_index, :] = (
    dolfin.interpolate(
        dolfin.Expression("x[0] > L / 2 ? 0.0 : g_Kr", g_Kr=g_Kr_value, L=L, degree=0),
        V_ode,
    )
    .vector()
    .get_local()
)

# and similar for the conductance $g_Ks$

g_Ks_index = model["parameter_index"]("g_Ks")
g_Ks_value = parameters[g_Ks_index]
parameters_ode[g_Ks_index, :] = (
    dolfin.interpolate(
        dolfin.Expression("x[0] > L / 2 ? 0.0 : g_Ks", g_Ks=g_Ks_value, L=L, degree=0),
        V_ode,
    )
    .vector()
    .get_local()
)

# Finally we set up the models

# +
pde = beat.MonodomainModel(
    time=time, mesh=mesh, M=D.magnitude, I_s=I_s, C_m=Cm.magnitude
)
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

fname = (outdir / "V.xdmf").as_posix()
beat.postprocess.delete_old_file(fname)


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


# -

# and solve it

t = 0.0
save_freq = int(1.0 / dt)
i = 0
done_stimulating = False
while t < end_time + 1e-12:
    # Make sure to save at the same time steps that is used by Ambit

    if i % save_freq == 0:
        save(t)
    if t > 1000 and not done_stimulating:
        ode.parameters[model["parameter_index"]("stim_amplitude"), :] = 0.0
        I_s.assign(0.0)
        done_stimulating = True

    solver.step((t, t + dt))
    i += 1
    t += dt


# Note that we also here stop stimulating after 1000 ms
#


def post_process(mesh, dx, outdir):
    V = dolfin.Function(dolfin.FunctionSpace(mesh, "CG", 1))
    xdmffile = (outdir / "V.xdmf").as_posix()

    cellnr = [0, 25, 50, 75, 100, 125, 150, 175, 200]
    cellnr = np.arange(0, 200, 10)
    times = beat.postprocess.load_timesteps_from_xdmf(xdmffile)
    t = np.array(list(times.values()))

    p1 = 25 * dx
    p2 = 175 * dx
    dp = 150 * dx * beat.units.ureg("cm")
    tp1 = np.inf
    tp2 = np.inf

    if not (outdir / "traces.npy").is_file():
        traces = np.zeros((len(times), len(cellnr)))
        with dolfin.XDMFFile(mesh.mpi_comm(), xdmffile) as xdmf:
            for i, ti in enumerate(times):
                xdmf.read_checkpoint(V, "V", ti)
                for j, cell in enumerate(cellnr):
                    traces[i, j] = V(cell * dx)

                if V(p1) > 0.0 and tp1 == np.inf:
                    tp1 = ti
                if V(p2) > 0.0 and tp2 == np.inf:
                    tp2 = ti

        np.save(outdir / "traces.npy", traces)
        np.save(outdir / "times.npy", t)
        (outdir / "cv.text").write_text(f"{tp1} {tp2}")

    traces = np.load(outdir / "traces.npy")
    t = np.load(outdir / "times.npy")
    tp1, tp2 = np.loadtxt(outdir / "cv.text")
    tp1 *= beat.units.ureg("ms")
    tp2 *= beat.units.ureg("ms")

    if not np.isclose(tp1, tp2):
        cv = dp / (tp2 - tp1)
        print(
            f"Conduction velocity:: {cv.to('cm/ms').magnitude} cm/ms "
            f"= {cv.to('m/s').magnitude} m/s"
        )

    # Plot 3D plot for all traces
    fig, ax = plt.subplots()
    for i, cell_index in enumerate(cellnr):
        color = "k" if cell_index < 100 else "m"
        ax.plot(t, cell_index / 3 * np.ones_like(t) + traces[:, i], color=color)
    ax.set_xlabel("Time (ms)")
    ax.set_yticks([-22, -55, -85])
    ax.set_yticklabels([200, 100, 1])
    fig.text(x=0.03, y=0.17, s="Cell number", rotation=90)

    ax.set_ylim(-90, 120)
    fig.savefig(outdir / "V_3d.png")


post_process(mesh, dx, outdir)


# ```{figure} ../docs/_static/pvc.png
# ---
# name: pvc
# ---
# Resulting voltage traces for different cells in the cable. Black cells represent normal cells, while colored cells represents cells with reduced conductance.
# ```

# ```{bibliography}
# ```
