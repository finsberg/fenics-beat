# # Niederer benchmark
# In this example we will use the same setup as in the Niederer benchmark
# > Niederer SA, Kerfoot E, Benson AP, Bernabeu MO, Bernus O, Bradley C,
#   Cherry EM, Clayton R, Fenton FH, Garny A, Heidenreich E, Land S, Maleckar M,
#   Pathmanathan P, Plank G, Rodríguez JF, Roy I, Sachse FB, Seemann G, Skavhaug O,
#   Smith NP. Verification of cardiac tissue electrophysiology simulators using an
#   N-version benchmark. Philos Trans A Math Phys Eng Sci. 2011 Nov 13;369(1954):4331-51.
#   doi: 10.1098/rsta.2011.0139. PMID: 21969679; PMCID: PMC3263775.
#

from pathlib import Path
import json

import dolfin
import numpy as np
import numpy.typing as npt
import pyvista
import pint
import gotranx
import matplotlib.pyplot as plt


import beat
import beat.viz

ureg = pint.UnitRegistry()


def setup_initial_conditions() -> npt.NDArray:
    ic = {
        "V": -85.23,  # mV
        "Xr1": 0.00621,
        "Xr2": 0.4712,
        "Xs": 0.0095,
        "m": 0.00172,
        "h": 0.7444,
        "j": 0.7045,
        "d": 3.373e-05,
        "f": 0.7888,
        "f2": 0.9755,
        "fCass": 0.9953,
        "s": 0.999998,
        "r": 2.42e-08,
        "Ca_i": 0.000126,  # millimolar
        "R_prime": 0.9073,
        "Ca_SR": 3.64,  # millimolar
        "Ca_ss": 0.00036,  # millimolar
        "Na_i": 8.604,  # millimolar
        "K_i": 136.89,  # millimolar
    }
    values = model.init_state_values(**ic)
    return values


def setup_geometry(Lx, Ly, Lz, dx):
    mesh = dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
        int(np.rint((Lz / dx))),
    )
    return mesh


def define_stimulus(mesh, chi, C_m, time, mesh_unit):
    S1_marker = 1
    # L = 1.5
    L = 1.5 * ureg("mm").to(mesh_unit).magnitude
    S1_subdomain = dolfin.CompiledSubDomain(
        "x[0] <= L + DOLFIN_EPS && x[1] <= L + DOLFIN_EPS && x[2] <= L + DOLFIN_EPS",
        L=L,
    )
    S1_markers = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    S1_subdomain.mark(S1_markers, S1_marker)

    # Define stimulation (NB: region of interest carried by the mesh
    # and assumptions in cbcbeat)
    duration = 2.0
    A = 50000.0 * ureg("uA / cm**3")

    amplitude = (A / (chi * C_m)).to("mV / ms").magnitude

    I_s = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=0.0,
        duration=duration,
        amplitude=amplitude,
        degree=0,
    )

    dx = dolfin.Measure("dx", domain=mesh, subdomain_data=S1_markers)(S1_marker)
    return beat.base_model.Stimulus(dz=dx, expr=I_s)


def define_conductivity_tensor(chi, C_m, mesh_unit):
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = 0.17 * ureg("S/m")
    sigma_it = 0.019 * ureg("S/m")
    sigma_el = 0.62 * ureg("S/m")
    sigma_et = 0.24 * ureg("S/m")

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(C_m * chi)
    s_l = (sigma_l / (C_m * chi)).to(f"{mesh_unit}**2/ms").magnitude
    s_t = (sigma_t / (C_m * chi)).to(f"{mesh_unit}**2/ms").magnitude

    # Define conductivity tensor
    M = dolfin.as_tensor(((s_l, 0, 0), (0, s_t, 0), (0, 0, s_t)))

    return M


# +
dx=0.5
dt=0.05
T=200.0

here = Path.cwd()

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

import tentusscher_panfilov_2006_epi_cell as model

fun = model.forward_generalized_rush_larsen
init_states = setup_initial_conditions()
parameters = model.init_parameter_values(stim_amplitude=0.0)

mesh_unit = "mm"

Lx = 20.0 * ureg("mm").to(mesh_unit).magnitude
Ly = 7 * ureg("mm").to(mesh_unit).magnitude
Lz = 3 * ureg("mm").to(mesh_unit).magnitude
dx_mm = dx * ureg("mm").to(mesh_unit).magnitude

mesh = setup_geometry(
    Lx=Lx,
    Ly=Ly,
    Lz=Lz,
    dx=dx_mm,
)
ode_space = dolfin.FunctionSpace(mesh, "Lagrange", 1)

pyvista.start_xvfb()
plotter = pyvista.Plotter()
topology, cell_types, x = beat.viz.create_vtk_structures(ode_space)
grid = pyvista.UnstructuredGrid(topology, cell_types, x)
plotter.add_mesh(grid, show_edges=True)

if not pyvista.OFF_SCREEN:
    plotter.show()
else:
    figure_as_array = plotter.screenshot("niederer-benchmark-geomeetry.png")

# +
# Surface to volume ratio
chi = 1400 * ureg("cm**-1")
# # Membrane capacitance
C_m = 1.0 * ureg("uF/cm**2")

time = dolfin.Constant(0.0)
I_s = define_stimulus(mesh=mesh, chi=chi, C_m=C_m, time=time, mesh_unit=mesh_unit)

M = define_conductivity_tensor(chi, C_m, mesh_unit=mesh_unit)

params = {"preconditioner": "sor", "use_custom_preconditioner": False}

pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
ode = beat.odesolver.DolfinODESolver(
    v_ode=dolfin.Function(ode_space),
    v_pde=pde.state,
    fun=fun,
    init_states=init_states,
    parameters=parameters,
    num_states=len(init_states),
    v_index=model.state_index("V"),
)

# +
solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode, theta=0.0)
output_dir = Path("output-niederer-benchmark")
output_dir.mkdir(exist_ok=True)
filename = output_dir / f"results-{dt}-{dx}.xdmf"
filename.unlink(missing_ok=True)
filename.with_suffix(".h5").unlink(missing_ok=True)

xdmf = dolfin.XDMFFile(mesh.mpi_comm(), filename.as_posix())
xdmf.parameters["functions_share_mesh"] = True
xdmf.parameters["rewrite_function_mesh"] = False
xdmf.write(mesh)

points = {
    "P1": (0, 0, 0),
    "P2": (0.0, Ly, 0.0),
    "P3": (Lx, 0.0, 0.0),
    "P4": (Lx, Ly, 0.0),
    "P5": (0.0, 0.0, Lz),
    "P6": (0.0, Ly, Lz),
    "P7": (Lx, 0.0, Lz),
    "P8": (Lx, Ly, Lz),
    "P9": (Lx / 2, Ly / 2, Lz / 2),
}

# +
activation_times = {p: -1.0 for p in points}
save_freq = int(1.0 / dt)
i = 0
plotter_voltage = pyvista.Plotter()
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
viridis = plt.get_cmap("viridis")
plotter_voltage.add_mesh(
    grid,
    show_edges=True,
    lighting=False,
    cmap=viridis,
    scalar_bar_args=sargs,
    clim=[-90.0, 40.0],
)
grid.point_data["V"] = solver.pde.state.vector().get_local()
gif_file = Path("niederer_benchmark.gif")
gif_file.unlink(missing_ok=True)
plotter_voltage.open_gif(gif_file.as_posix())

T = 10
# T = 100  # Change to 100 to reproduce Niederer benchmark
t = 0.0
while t < T + 1e-12 and any(at < 0.0 for at in activation_times.values()):
    v = solver.pde.state.vector().get_local()
    if i % save_freq == 0:
        print(f"Solve for {t=:.2f}, {v.max() =}, {v.min() = }")
        print(activation_times)
        xdmf.write_checkpoint(
            solver.pde.state,
            "V",
            float(t),
            dolfin.XDMFFile.Encoding.HDF5,
            True,
        )
        grid.point_data["V"] = solver.pde.state.vector().get_local()
        plotter_voltage.write_frame()
    solver.step((t, t + dt))

    for p in points:
        value = beat.utils.peval(solver.pde.state, points[p])
        if value > 0.0 and activation_times[p] < 0.0:
            activation_times[p] = t
    i += 1
    t += dt

mesh.mpi_comm().Barrier()
plotter_voltage.close()
# -

# ![niederer](niederer_benchmark.gif "niederer")

if mesh.mpi_comm().rank == 0:
    # Save activation times
    activation_times["dx"] = dx
    activation_times["dt"] = dt
    at_file_name = output_dir / "activation_times.json"
    if at_file_name.is_file():
        all_at = json.loads(at_file_name.read_text())
    else:
        all_at = []
    all_at.append(activation_times)
    at_file_name.write_text(json.dumps(all_at, indent=2))

# The activation times are saved in the file `output-niederer-benchmark/activation_times.json`.
# The file contains a list of dictionaries, each dictionary contains the activation times for a specific dx and dt.
# Here are the activation times for the different dx and dt:
#
# |    |   dx |    dt |    P1 |     P2 |     P3 |     P4 |     P5 |     P6 |     P7 |     P8 |     P9 |
# |---:|-----:|------:|------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|-------:|
# |  0 |  0.5 | 0.05  | 1.25  | 51.1   | 34.9   | 58.9   | 14.1   | 49.5   | 34     | 56.65  | 26.05  |
# |  1 |  0.5 | 0.01  | 1.22  | 50.85  | 33.96  | 58.05  | 13.98  | 49.36  | 33.07  | 55.91  | 25.64  |
# |  2 |  0.5 | 0.005 | 1.215 | 50.775 | 33.825 | 57.96  | 13.97  | 49.345 | 32.945 | 55.825 | 25.595 |
# |  3 |  0.2 | 0.05  | 1.25  | 29.7   | 32.9   | 40.2   |  9.55  | 30     | 32.95  | 39.9   | 18.9   |
# |  4 |  0.2 | 0.01  | 1.24  | 29.09  | 31.25  | 38.66  |  9.34  | 29.4   | 31.29  | 38.42  | 18.14  |
# |  5 |  0.2 | 0.005 | 1.235 | 29.015 | 31.05  | 38.475 |  9.315 | 29.32  | 31.08  | 38.235 | 18.045 |
# |  6 |  0.1 | 0.05  | 1.25  | 26.85  | 33.3   | 40.35  |  8.4   | 27.5   | 33.85  | 40.55  | 18.95  |
# |  7 |  0.1 | 0.01  | 1.23  | 25.64  | 31.46  | 38.08  |  8.03  | 26.24  | 31.94  | 38.21  | 17.95  |
# |  8 |  0.1 | 0.005 | 1.225 | 25.5   | 31.26  | 37.81  |  7.99  | 26.09  | 31.72  | 37.93  | 17.835 |
