# # Stimulation with Purkinje system
# In this demo we show how to simulate endocardial stimulation using a Purkinje system for activation.
#
# We use the package [fractal_tree](https://gitub.com/finsberg/fractal_tree) to generate the Purkinje system. First let import the neccessary packages

from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import dolfin
import pyvista
import meshio
import cardiac_geometries
import beat
import beat.viz
import gotranx
from fractal_tree import generate_fractal_tree, FractalTreeParameters, Mesh


# We create a function that generates the BiV geometry that we want to use for this demo.


def get_data(datadir="data_endocardial_stimulation"):
    datadir = Path(datadir)
    msh_file = datadir / "biv_ellipsoid.msh"
    if not msh_file.is_file():
        cardiac_geometries.create_biv_ellipsoid(
            datadir,
            char_length=0.15,  # Reduce this value to get a finer mesh
            center_lv_y=0.2,
            center_lv_z=0.0,
            a_endo_lv=5.0,
            b_endo_lv=2.2,
            c_endo_lv=2.2,
            a_epi_lv=6.0,
            b_epi_lv=3.0,
            c_epi_lv=3.0,
            center_rv_y=1.0,
            center_rv_z=0.0,
            a_endo_rv=6.0,
            b_endo_rv=2.5,
            c_endo_rv=2.7,
            a_epi_rv=8.0,
            b_epi_rv=5.5,
            c_epi_rv=4.0,
            create_fibers=True,
        )

    return cardiac_geometries.geometry.Geometry.from_folder(datadir)


# We have a function that can load the solution from a file.


def load_from_file(heart_mesh, xdmffile, key="v", stop_index=None):
    V = dolfin.FunctionSpace(heart_mesh, "Lagrange", 1)
    v = dolfin.Function(V)

    timesteps = beat.postprocess.load_timesteps_from_xdmf(xdmffile)
    with dolfin.XDMFFile(Path(xdmffile).as_posix()) as f:
        for i, t in timesteps.items():
            f.read_checkpoint(v, key, i)
            yield v.copy(deepcopy=True), t


# and we have a function to compute the pseudo ECG


def compute_ecg_recovery():
    datadir = Path("data_endocardial_stimulation")
    xdmffile = datadir / "state.xdmf"
    data = get_data(datadir=datadir)

    # https://litfl.com/ecg-lead-positioning/
    vs = load_from_file(data.mesh, xdmffile, key="V")

    leads = dict(
        RA=(-15.0, 0.0, -10.0),
        LA=(4.0, -12.0, -7.0),
        RL=(0.0, 20.0, 3.0),
        LL=(17.0, 11.0, 7.0),
        V1=(-3.0, 4.0, -9.0),
        V2=(0.0, 2.0, -8.0),
        V3=(3.0, 1.0, -8.0),
        V4=(6.0, 1.0, -6.0),
        V5=(10.0, 2.0, 0.0),
        V6=(10.0, -6.0, 2.0),
    )

    fname = datadir / "extracellular_potential.npy"
    if 1:  # not fname.is_file():
        phie = defaultdict(list)
        time = []
        for v, t in vs:
            time.append(t)
            for name, point in leads.items():
                phie[name].append(
                    beat.ecg.ecg_recovery(v=v, mesh=data.mesh, sigma_b=1.0, point=point)
                )
        np.save(fname, {"phie": phie, "time": time})

    phie_time = np.load(fname, allow_pickle=True).item()
    phie = phie_time["phie"]
    time = phie_time["time"]

    fig, ax = plt.subplots(2, 5, sharex=True, figsize=(12, 8))
    for i, (name, values) in enumerate(phie.items()):
        axi = ax.flatten()[i]
        axi.plot(time, values)
        axi.set_title(name)
    fig.savefig(datadir / "extracellular_potential.png")

    ecg = beat.ecg.Leads12(**{k: np.array(v) for k, v in phie.items()})
    fig, ax = plt.subplots(3, 4, sharex=True, figsize=(12, 8))
    for i, name in enumerate(
        [
            "I",
            "II",
            "III",
            "aVR",
            "aVL",
            "aVF",
            "V1_",
            "V2_",
            "V3_",
            "V4_",
            "V5_",
            "V6_",
        ]
    ):
        y = getattr(ecg, name)
        axi = ax.flatten()[i]
        axi.plot(time, y)
        axi.set_title(name)
    fig.savefig(datadir / "ecg_12_leads.png")


# First let use use [gotranx](https://finsberg.github.io/gotranx/) to generate code for the ORd model.

# +
here = Path.cwd()

model_path = here / "ORdmm_Land.py"
if not model_path.is_file():
    ode = gotranx.load_ode(here / ".." / "odes" / "ORdmm_Land.ode")
    code = gotranx.cli.gotran2py.get_code(
        ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
    )
    model_path.write_text(code)

import ORdmm_Land

model = ORdmm_Land.__dict__
# -

# Next we generate the geometry that we will use for this demo and save it to the folder `data_endocardial_stimulation`.

datadir = here / "data_endocardial_stimulation"
data = get_data(datadir=datadir)

# We can now generate the Purkinje system for the LV. We start by loading the mesh that is generated by `GMSH` and extracting the endocardium of the LV. We first read the gmsh file

msh_file = datadir / "biv_ellipsoid.msh"

# We load file with meshio

msh = meshio.read(msh_file)

# and we load the geometry using `cardiac-geometries`

geo = cardiac_geometries.geometry.Geometry.from_folder(datadir)

# We can take a short look at the geometry in pyvista

V = dolfin.FunctionSpace(geo.mesh, "Lagrange", 1)
# pyvista.start_xvfb()
# plotter_markers = pyvista.Plotter()
# topology, cell_types, x = beat.viz.create_vtk_structures(V)
# grid = pyvista.UnstructuredGrid(topology, cell_types, x)
# plotter_markers.add_mesh(grid, show_edges=True)
# plotter_markers.view_zy()
# plotter_markers.camera.zoom(4)
# if not pyvista.OFF_SCREEN:
#     plotter_markers.show()
# else:
#     figure = plotter_markers.screenshot("biv_mesh_purkinje.png")


# Next we extract the coordinates of the mesh

verts = geo.mesh.coordinates()

# and then we extract the cells for the LV endocardium

tag = msh.field_data["ENDO_LV"][0]
inds = [i for i, x in enumerate(msh.cell_data["gmsh:physical"]) if x[0] == tag]
connectivity = np.vstack([msh.cells[i].data for i in inds])

# To generate the purkinje tree we must first define the initial node for the purkinje system in the LV. Here we have uses paraview to an approximate coordinate for this which is located near the base on the endocardium.

init_node = [0, 2.4, 0.19]

# Now that we have the approximate coordinate, we need to find the closest node in the mesh. We can do this by selecting the point with the shortest distance to this node.

index = np.linalg.norm(np.subtract(verts, init_node), axis=1).argmin()

# We create the mesh, and pass in the initial node
lv_mesh = Mesh(verts=verts, connectivity=connectivity, init_node=verts[index, :])

# We also set the number of generations to 15 and we set the initial direction
# to point in the positive $x$-direction which is the direction from the base
# towards the apex (this can also be found by visualizing the mesh in Paraview).
# We also specify the initial length to be 3. This is about the length of the
# septal wall, so that the first branch will represent the His bundle.
# We set the length to 0.25 which will be the length each successive branch.
#
# First we set a seed so that we know that the results are reproducible

np.random.seed(1)

# Assume conduction velocty of 3 m/s = 3 mm/ms

# +
conduction_velocity = 3.0
duration = 2.0
start = 1.0
dofs = V.tabulate_dof_coordinates()

param = FractalTreeParameters(
    filename=datadir / "lv_tree",
    init_length=7.0,
    N_it=15,
    length=0.5,
    initial_direction=np.array([1, 0, 0]),
)
# -

# Next we create the Purkinje networks for the LV
# +
lv_tree = generate_fractal_tree(lv_mesh, param)

# pyvista.start_xvfb()
# lv_tree_vtu = pyvista.read((datadir / "lv_tree").with_suffix(".vtu"))
# plotter = pyvista.Plotter()
# plotter.add_mesh(lv_tree_vtu, show_edges=True)

# if not pyvista.OFF_SCREEN:
#     plotter.show()
# else:
#     figure = plotter.screenshot("lv_tree.png")
# # -


# We can now compute the distance from the root node to each node in the tree

dist_lv = np.zeros(lv_tree.nodes.shape[0])
for line in lv_tree.lines:
    dist_lv[line[1]] += dist_lv[line[0]] + np.linalg.norm(
        lv_tree.nodes[line[1], :] - lv_tree.nodes[line[0], :]
    )

# Then we extract the end nodes of the tree and compute the distance from the root

lv_end_nodes = np.array(lv_tree.end_nodes, dtype=int)
lv_end_nodes_pos = lv_tree.nodes[lv_end_nodes]
lv_end_nodes_dist = dist_lv[lv_end_nodes]

# We can now compute the delay for each end node given the conduction velocity

delay_lv = lv_end_nodes_dist / conduction_velocity

# Finally, we create the expression for the stimulation protocol

delay_expr_lv = []
for pos, d in zip(lv_end_nodes_pos, delay_lv):
    closest_dof = np.linalg.norm(dofs - pos, axis=1).argmin()
    delay_expr_lv.append(
        f"(time > {start + d} && time < {start + d + duration} && near(x[0], {dofs[closest_dof, 0]:.2f}, tol) && near(x[1], {dofs[closest_dof, 1]:.2f}, tol) && near(x[2], {dofs[closest_dof, 2]:.2f}, tol))"
    )

# Repeat for RV

tag = msh.field_data["ENDO_RV"][0]
inds = [i for i, x in enumerate(msh.cell_data["gmsh:physical"]) if x[0] == tag]
connectivity = np.vstack([msh.cells[i].data for i in inds])

# and defining the initial node (again, we use Paraview to find this location).

init_node = [0, 3.19, 0.1]
index = np.linalg.norm(np.subtract(verts, init_node), axis=1).argmin()
rv_mesh = Mesh(verts=verts, connectivity=connectivity, init_node=verts[index, :])

# Set a new seed (note that if you don't like the generated system you can try a different seed)
# +
np.random.seed(1234)

param = FractalTreeParameters(
    filename=datadir / "rv_tree",
    init_length=9.7,
    N_it=15,
    length=0.5,
    initial_direction=np.array([1, 0, 0]),
)

rv_tree = generate_fractal_tree(rv_mesh, param)
# -

# We can now compute the distance from the root node to each node in the tree

# +
dist_rv = np.zeros(rv_tree.nodes.shape[0])
for line in rv_tree.lines:
    dist_rv[line[1]] += dist_rv[line[0]] + np.linalg.norm(
        rv_tree.nodes[line[1], :] - rv_tree.nodes[line[0], :]
    )

rv_end_nodes = np.array(rv_tree.end_nodes, dtype=int)
rv_end_nodes_pos = rv_tree.nodes[rv_end_nodes]
rv_end_nodes_dist = dist_rv[rv_end_nodes]

delay_rv = rv_end_nodes_dist / conduction_velocity

delay_expr_rv = []
for pos, d in zip(rv_end_nodes_pos, delay_rv):
    closest_dof = np.linalg.norm(dofs - pos, axis=1).argmin()
    delay_expr_rv.append(
        f"(time > {start + d} && time < {start + d + duration} && near(x[0], {dofs[closest_dof, 0]:.2f}, tol) && near(x[1], {dofs[closest_dof, 1]:.2f}, tol) && near(x[2], {dofs[closest_dof, 2]:.2f}, tol))"
    )

expr = (
    f"time < {start + max(delay_lv.max(), delay_rv.max()) + duration} && ("
    + " || ".join(delay_expr_lv + delay_expr_rv)
    + ") ? stim_amp : 0.0"
)


# +
mesh_unit = "mm"

time = dolfin.Constant(0.0)
stim_expr = dolfin.Expression(
    expr,
    time=time,
    tol=0.05,
    stim_amp=0.01,
    degree=1,
)

markers = beat.utils.expand_layer_biv(
    V=V,
    mfun=data.ffun,
    endo_lv_marker=data.markers["ENDO_LV"][0],
    endo_rv_marker=data.markers["ENDO_RV"][0],
    epi_marker=data.markers["EPI"][0],
    endo_size=0.3,
    epi_size=0.3,
)


# +
init_states = {
    0: model["init_state_values"](),
    1: model["init_state_values"](),
    2: model["init_state_values"](),
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

# Surface to volume ratio
chi = 140.0 * beat.units.ureg("mm**-1")
# Membrane capacitance
C_m = 0.01 * beat.units.ureg("uF/mm**2")

I_s = beat.stimulation.Stimulus(dolfin.dP(domain=data.mesh), expr=stim_expr)


M = beat.conductivities.define_conductivity_tensor(
    chi=chi, f0=data.f0, g_il=0.17, g_it=0.019, g_el=0.62, g_et=0.24
)

params = {"preconditioner": "sor", "use_custom_preconditioner": False}
pde = beat.MonodomainModel(
    time=time,
    C_m=C_m.to(f"uF/{mesh_unit}**2").magnitude,
    mesh=data.mesh,
    M=M,
    I_s=I_s,
    params=params,
)

ode = beat.odesolver.DolfinMultiODESolver(
    v_ode=dolfin.Function(V),
    v_pde=pde.state,
    markers=markers,
    num_states={i: len(s) for i, s in init_states.items()},
    fun=fun,
    init_states=init_states,
    parameters=parameters,
    v_index=v_index,
)

# # +
T = 3
# Change to 500 to simulate the full cardiac cycle
t = 0.0
dt = 0.05
solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

fname = (datadir / "state.xdmf").as_posix()
i = 0
while t < T + 1e-12:
    if i % 20 == 0:
        v = solver.pde.state.vector().get_local()
        print(f"Solve for {t=:.2f}, {v.max() =}, {v.min() = }")
        with dolfin.XDMFFile(dolfin.MPI.comm_world, fname) as xdmf:
            xdmf.parameters["functions_share_mesh"] = True
            xdmf.parameters["rewrite_function_mesh"] = False
            xdmf.write_checkpoint(
                solver.pde.state,
                "V",
                float(t),
                dolfin.XDMFFile.Encoding.HDF5,
                True,
            )
    solver.step((t, t + dt))
    i += 1
    t += dt
# # -

compute_ecg_recovery()
