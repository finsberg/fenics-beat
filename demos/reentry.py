# # Conduction velocity and ECG for slabs
# In this demo we will show how to compute conduction velocity and ECG for a Slab geometry.
#

from pathlib import Path
import beat.single_cell
import dolfin
import numpy as np

try:
    import ufl_legacy as ufl
except ImportError:
    import ufl

import gotranx
import beat
import beat.viz


def define_conductivity_tensor(
    chi,
    f0,
    V_dg,
    indices,
    g_il=0.17,
    g_it=0.019,
    g_el=0.62,
    g_et=0.24,
):
    s_l, s_t = beat.conductivities.get_harmonic_mean_conductivity(
        chi=chi, g_il=g_il, g_it=g_it, g_el=g_el, g_et=g_et
    )

    s_l_func = dolfin.Function(V_dg)
    s_l_func.vector()[:] = s_l
    s_l_func.vector()[indices] = 1e-7
    s_t_func = dolfin.Function(V_dg)
    s_t_func.vector()[:] = s_t
    s_t_func.vector()[indices] = 1e-7

    return s_l * ufl.outer(f0, f0) + s_t * (ufl.Identity(2) - ufl.outer(f0, f0))


def setup_geometry(dx, Lx, Ly):

    mesh = dolfin.RectangleMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0),
        dolfin.Point(Lx, Ly),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
    )

    return mesh


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

dx = 0.01 * beat.units.ureg("cm").to(mesh_unit).magnitude
L = 1.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
mesh = setup_geometry(Lx=L, Ly=L, dx=dx)


g_il = 0.16069
g_el = 0.625
g_it = 0.04258
g_et = 0.236

f0, s0 = get_microstructure(dimension, transverse)


data = beat.Geometry(
    mesh=mesh,
    f0=f0,
)
save_freq = round(save_every_ms / dt)


print("Running model")
# # Load the model
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

subdomain_data = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
subdomain_data.set_all(0)
marker = 1
dolfin.CompiledSubDomain(
    "std::pow((x[0] - 0.2), 2) + std::pow((x[1] - 0.5), 2) < 0.05", dx=dx
).mark(subdomain_data, 1)

I_s = beat.stimulation.define_stimulus(
    chi=chi,
    mesh=mesh,
    mesh_unit=mesh_unit,
    time=time,
    start=5.0,
    duration=5.0,
    amplitude=1.0,
    PCL=150.0,
    subdomain_data=subdomain_data,
    marker=marker,
)

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

fname = (results_folder / "V.xdmf").as_posix()
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

i = 0
while t < end_time + 1e-12:
    # Make sure to save at the same time steps that is used by Ambit

    if i % save_freq == 0:
        save(t)

    solver.step((t, t + dt))
    i += 1
    t += dt
