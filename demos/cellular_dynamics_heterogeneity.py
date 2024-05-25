# # Cellular dynamics heterogeneity
#
#

from typing import NamedTuple
from pathlib import Path
import beat.single_cell
import dolfin
import numpy as np

import gotranx
import beat
import pyvista
import beat.viz


class Geometry(NamedTuple):
    mesh: dolfin.Mesh
    ffun: dolfin.MeshFunction
    markers: dict[str, tuple[int, int]]
    f0: dolfin.Constant
    s0: dolfin.Constant


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
    data: Geometry,
    markers: dolfin.Function,
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
    # Load the model
    model_path = Path("mahajan.py")

    if not model_path.is_file():
        print("Generate code for cell model")
        here = Path.cwd()
        ode = gotranx.load_ode(
            here
            / ".."
            / "odes"
            / "mahajan_shiferaw_sato_baher_olcese_xie_yang_chen_restrepo_karma_garfinkel_qu_weiss_2008.ode"
        )
        code = gotranx.cli.gotran2py.get_code(
            ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
        )
        model_path.write_text(code)

    import mahajan

    model = mahajan.__dict__

    # Surface to volume ratio
    chi = 1400.0 * beat.units.ureg("cm**-1")

    # Membrane capacitance
    C_m = 1.0 * beat.units.ureg("uF/cm**2")

    with dolfin.XDMFFile((resultsdir / "markers.xdmf").as_posix()) as xdmf:
        xdmf.write(markers)

    fun = model["forward_generalized_rush_larsen"]

    y = model["init_state_values"]()

    time = dolfin.Constant(0.0)
    parameters = model["init_parameter_values"](stim_amplitude=0.0)

    I_s_expr = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=5.0,
        duration=2.0,
        amplitude=1.0,
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
    g_Kr_index = model["parameter_index"]("gkr")
    g_Kr_value = parameters[g_Kr_index]
    g_Kr_func = dolfin.interpolate(markers, V_ode)
    g_Kr_func.vector()[:] *= g_Kr_value
    parameters_ode[g_Kr_index, :] = g_Kr_func.vector().get_local()

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
    end_time = 300.0
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


results_folder = Path("results-cellular-dynamics-heterogeneity")
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

dx = 0.05 * beat.units.ureg("cm").to(mesh_unit).magnitude
L = 1.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
mesh = setup_geometry(Lx=L, Ly=L, dx=dx)

ffun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
ffun.set_all(0)
stim_domain = dolfin.CompiledSubDomain("x[0] <= DOLFIN_EPS")
marker = 1
stim_domain.mark(ffun, marker)

cfun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
cfun.set_all(0)
imp_region_0 = dolfin.CompiledSubDomain("x[1] < L / 4 + DOLFIN_EPS", L=L)
imp_region_1 = dolfin.CompiledSubDomain(
    "x[1] > L / 4 - DOLFIN_EPS && x[1] < L / 2 + DOLFIN_EPS", L=L
)
imp_region_2 = dolfin.CompiledSubDomain(
    "x[1] > L / 2 - DOLFIN_EPS && x[1] < 3 * L / 4 + DOLFIN_EPS", L=L
)
imp_region_3 = dolfin.CompiledSubDomain("x[1] > 3 * L / 4 - DOLFIN_EPS", L=L)

imp_region_0.mark(cfun, 1)
imp_region_1.mark(cfun, 2)
imp_region_2.mark(cfun, 3)
imp_region_3.mark(cfun, 4)

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

f0, s0 = get_microstructure(transverse)

markers = {"ENDO": (marker, 2)}

data = Geometry(
    mesh=mesh,
    ffun=ffun,
    markers=markers,
    f0=f0,
    s0=s0,
)
save_freq = round(save_every_ms / dt)
run_model(
    data=data,
    markers=cfun_func,
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
