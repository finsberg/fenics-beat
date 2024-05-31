# # Conduction velocity and ECG for slabs
# In this demo we will show how to compute conduction velocity and ECG for a Slab geometry.
#

from pathlib import Path
import beat.conductivities
import beat.postprocess
import beat.single_cell
import beat.stimulation
import dolfin
import matplotlib.pyplot as plt
import numpy as np
import gotranx
import beat
import pyvista
import beat.viz
from beat.geometry import Geometry


results_folder = Path("results-slab")
save_every_ms = 1.0
dimension = 2
transverse = False
end_time = 20.0
dt = 0.05
overwrite = False
stim_amp = 5000.0
mesh_unit = "cm"
dx = 0.05 * beat.units.ureg("cm").to(mesh_unit).magnitude
L = 1.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
mesh = beat.geometry.get_3D_slab_mesh(Lx=L, Ly=dx, Lz=dx, dx=dx / 5, dim=dimension)

ffun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
ffun.set_all(0)
stim_domain = dolfin.CompiledSubDomain("x[0] <= DOLFIN_EPS")
marker = 1
stim_domain.mark(ffun, marker)

cfun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
cfun.set_all(0)
endo = dolfin.CompiledSubDomain("x[0] < L / 3", L=L)
epi = dolfin.CompiledSubDomain("x[0] > 2 * L / 3", L=L)
endo.mark(cfun, 1)
epi.mark(cfun, 2)

V = dolfin.FunctionSpace(mesh, "CG", 1)


pyvista.start_xvfb()
plotter_markers = pyvista.Plotter()
topology, cell_types, x = beat.viz.create_vtk_structures(V)
grid = pyvista.UnstructuredGrid(topology, cell_types, x)
grid["markers"] = cfun.array()
plotter_markers.add_mesh(grid, show_edges=True)
if mesh.geometric_dimension() == 2:
    plotter_markers.view_xy()

if not pyvista.OFF_SCREEN:
    plotter_markers.show()
else:
    figure_as_array = plotter_markers.screenshot("markers.png")

#
# Interpolate meshfunction to a CG 1 function
#

cfun_DG = dolfin.Function(dolfin.FunctionSpace(mesh, "DG", 0))
cfun_DG.vector()[:] = cfun.array()
cfun_func = dolfin.Function(V)
cfun_func.interpolate(cfun_DG)

g_il = 0.16069
g_el = 0.625
g_it = 0.04258
g_et = 0.236

f0, s0, n0 = beat.geometry.get_3D_slab_microstructure(transverse)

markers = {"ENDO": (marker, 2)}

data = Geometry(
    mesh=mesh,
    ffun=ffun,
    markers=markers,
    f0=f0,
    s0=s0,
    n0=n0,
)
save_freq = round(save_every_ms / dt)


print("Running model")
# Load the model
model_path = Path("ORdmm_Land.py")
if not model_path.is_file():
    print("Generate code for cell model")
    here = Path.cwd()
    ode = gotranx.load_ode(here / ".." / "odes" / "ORdmm_Land.ode")
    code = gotranx.cli.gotran2py.get_code(
        ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
    )
    model_path.write_text(code)

import ORdmm_Land

model = ORdmm_Land.__dict__

# Surface to volume ratio

chi = 1400.0 * beat.units.ureg("cm**-1")

# Membrane capacitance

C_m = 1.0 * beat.units.ureg("uF/cm**2")


print("Get steady states")
nbeats = 2  # Should be set to at least 200
init_states = {
    0: beat.single_cell.get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=2),
        outdir=results_folder / "mid",
        BCL=1000,
        nbeats=nbeats,
        track_indices=[model["state_index"]("v"), model["state_index"]("cai")],
        dt=0.05,
    ),
    1: beat.single_cell.get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=0),
        outdir=results_folder / "endo",
        BCL=1000,
        nbeats=nbeats,
        track_indices=[
            model["state_index"]("v"),
            model["state_index"]("cai"),
            model["state_index"]("nai"),
        ],
        dt=0.05,
    ),
    2: beat.single_cell.get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=1),
        outdir=results_folder / "epi",
        BCL=1000,
        nbeats=nbeats,
        track_indices=[model["state_index"]("v"), model["state_index"]("cai")],
        dt=0.05,
    ),
}
# endo = 0, epi = 1, M = 2
parameters = {
    0: model["init_parameter_values"](amp=0.0, celltype=2),
    1: model["init_parameter_values"](amp=0.0, celltype=0),
    2: model["init_parameter_values"](amp=0.0, celltype=1),
}
fun = {
    0: model["forward_generalized_rush_larsen"],
    1: model["forward_generalized_rush_larsen"],
    2: model["forward_generalized_rush_larsen"],
}
v_index = {
    0: model["state_index"]("v"),
    1: model["state_index"]("v"),
    2: model["state_index"]("v"),
}

time = dolfin.Constant(0.0)
I_s = beat.stimulation.define_stimulus(
    mesh=data.mesh,
    chi=chi,
    time=time,
    subdomain_data=data.ffun,
    marker=markers["ENDO"][0],
    amplitude=stim_amp,
    mesh_unit=mesh_unit,
)

M = beat.conductivities.define_conductivity_tensor(
    chi,
    f0=data.f0,
    g_il=g_il,
    g_it=g_it,
    g_el=g_el,
    g_et=g_et,
)

params = {"preconditioner": "sor", "use_custom_preconditioner": False}
pde = beat.MonodomainModel(
    time=time,
    mesh=data.mesh,
    M=M,
    I_s=I_s,
    params=params,
    C_m=C_m.to(f"uF/{mesh_unit}**2").magnitude,
)

V_ode = dolfin.FunctionSpace(data.mesh, "Lagrange", 1)
ode = beat.odesolver.DolfinMultiODESolver(
    v_ode=dolfin.Function(V_ode),
    v_pde=pde.state,
    markers=cfun_func,
    num_states={i: len(s) for i, s in init_states.items()},
    fun=fun,
    init_states=init_states,
    parameters=parameters,
    v_index=v_index,
)
t = 0.0
solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

fname = (results_folder / "V.xdmf").as_posix()
beat.postprocess.delete_old_file(fname)


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


i = 0
while t < end_time + 1e-12:
    # Make sure to save at the same time steps that is used by Ambit

    if i % save_freq == 0:
        save(t)

    solver.step((t, t + dt))
    i += 1
    t += dt


# ## Compute conduction velocity

threshold = 0.0
x0 = L * 0.25
x1 = L * 0.75
if mesh.geometry().dim() == 2:
    p1 = (x0, dx * 0.5)
    p2 = (x1, dx * 0.5)
else:
    p1 = (x0, dx * 0.5, dx * 0.5)  # type: ignore
    p2 = (x1, dx * 0.5, dx * 0.5)  # type: ignore

V = dolfin.FunctionSpace(mesh, "CG", 1)
v = dolfin.Function(V)


plotter_voltage = pyvista.Plotter()
grid = pyvista.UnstructuredGrid(topology, cell_types, x)
grid.point_data["V"] = v.vector().get_local()
viridis = plt.get_cmap("viridis")
sargs = dict(
    title_font_size=25,
    label_font_size=20,
    fmt="%.2e",
    color="black",
    position_x=0.1,
    position_y=0.8,
    width=0.8,
    height=0.1,
)

# +
plotter_voltage.add_mesh(
    grid,
    show_edges=True,
    lighting=False,
    cmap=viridis,
    scalar_bar_args=sargs,
    clim=[-90.0, 40.0],
)
fname = (results_folder / "V.xdmf").as_posix()
times = beat.postprocess.load_timesteps_from_xdmf(fname)
t1 = np.inf
t2 = np.inf
gif_file = Path("voltage_slab_time.gif")
gif_file.unlink(missing_ok=True)
plotter_voltage.open_gif(gif_file.as_posix())
for i, t in times.items():
    with dolfin.XDMFFile(mesh.mpi_comm(), fname) as xdmf:
        xdmf.read_checkpoint(v, "V", i)
        print(f"Read {t=:.2f}, {v(p1) =}, {v(p2) = }")

    grid.point_data["V"] = v.vector().get_local()
    plotter_voltage.write_frame()

    if v(p1) > threshold:
        t1 = min(t, t1)
    if v(p2) > threshold:
        t2 = min(t, t2)

plotter_voltage.close()
# -

# ![volt](voltage_slab_time.gif "volt")

if not np.isclose(t1, t2):
    cv = (x1 - x0) / (t2 - t1) * beat.units.ureg(f"{mesh_unit}/ms")
    msg = (
        f"Conduction velocity = {cv.magnitude:.3f} mm/ms or "  #
        f" {cv.to('m/s').magnitude:.3f} m/s or "  #
        f" {cv.to('cm/s').magnitude:.3f} cm/s"  #
    )
    print(msg)

x0 = L * 2.0
if mesh.geometry().dim() == 2:
    p = (x0, dx * 0.5)
else:
    p = (x0, dx * 0.5, dx * 0.5)  # type: ignore

V = dolfin.FunctionSpace(mesh, "CG", 1)
v = dolfin.Function(V)

phie = []
for i, t in times.items():
    with dolfin.XDMFFile(mesh.mpi_comm(), fname) as xdmf:
        xdmf.read_checkpoint(v, "V", i)
        print(f"Read {t=:.2f}")
        phie.append(beat.ecg.ecg_recovery(v=v, mesh=mesh, sigma_b=1.0, point=p))

fig, ax = plt.subplots(1, 1, figsize=(12, 8))
ax.plot(times.values(), phie)
ax.set_title("ECG recovery")
fig.savefig(results_folder / "ecg_recovery.png")
