# # Fibrotic region
# In this demo we will show how to model a fibrotic region in a 2D sheet of cardiac tissue. This demo is inspired by the demo in https://opencarp.org/documentation/examples/02_ep_tissue/21_reentry_induction

# First we do the necessary imports

# +
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


# -


# We will define a function to define the conductivity tensor which will
# modify the conductivity tensor in the fibrotic region of the tissue and basically the the conductivities to $10^{-7}$ in these elements.
#


def define_conductivity_tensor(
    chi,
    f0,
    V_dg,
    indices,
    g_il,
    g_it,
    g_el,
    g_et,
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

    return s_l_func * ufl.outer(f0, f0) + s_t_func * (
        ufl.Identity(2) - ufl.outer(f0, f0)
    )


# Now just define some general parameters for the simulation

results_folder = Path("results-fibrotic-region")
save_every_ms = 1.0
dimension = 2
transverse = False
# Increase this to see more interesting dynamics
end_time = 2.0
dt = 0.05
save_freq = round(save_every_ms / dt)
overwrite = False
stim_amp = 50_000.0
mesh_unit = "cm"
dx = 0.4 * beat.units.ureg("mm").to(mesh_unit).magnitude
L = 5.0 * beat.units.ureg("cm").to(mesh_unit).magnitude
data = beat.geometry.get_2D_slab_geometry(
    Lx=L,
    Ly=L,
    dx=dx,
    transverse=transverse,
)


# +
print("Running model")
# # Load the model
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

# Define model parameters for conductivities
conductivities = beat.conductivities.default_conductivities("Niederer")
# and the membrane capacitance
C_m = 1.0 * beat.units.ureg("uF/cm**2")

# Gather the defails initial states and parameters and set the stimulus amplitude to zero as this will be provided at the PDE level
fun = model["forward_generalized_rush_larsen"]
y = model["init_state_values"]()
time = dolfin.Constant(0.0)
parameters = model["init_parameter_values"](stim_amplitude=0.0)
# Next we define the point where we want to stimulate and choose a small region close to the left boundary

subdomain_data = dolfin.MeshFunction("size_t", data.mesh, data.mesh.topology().dim())
subdomain_data.set_all(0)
marker = 1
dolfin.CompiledSubDomain("x[0] < 2*dx", dx=dx).mark(subdomain_data, 1)
dolfin.File("subdomain_data.pvd") << subdomain_data

# and define the stimulus to stimulate every 350 ms with a duration of 5 ms

I_s = beat.stimulation.define_stimulus(
    chi=conductivities["chi"],
    mesh=data.mesh,
    mesh_unit=mesh_unit,
    time=time,
    start=5.0,
    duration=5.0,
    amplitude=stim_amp,
    PCL=350.0,
    subdomain_data=subdomain_data,
    marker=marker,
)

# Next we will define the fibrotic region is defined as a region in the center of radius 1.42 cm, i.e
#
# ```{math}
# (x - L/2)^2 + (y - L/2)^2 < 1^2
# ```
#

fibrosis = dolfin.MeshFunction("size_t", data.mesh, data.mesh.topology().dim())
fibrosis.set_all(0)
marker = 1
dolfin.CompiledSubDomain(
    "std::pow((x[0] - L/2),2) + std::pow((x[1] - L/2), 2) < std::pow(1.42, 2)", L=L
).mark(fibrosis, marker)

# Inside this region we will have 70% of the tissue with reduced conductance for the `g_K1`, `g_Na` and `g_Ca_L` and the remaining 30% will have a conductivity of $10^{-7}$.
#
# We first find the indices and make 70/30 split

indices = np.where(fibrosis.array() == 1)[0]
np.random.shuffle(indices)
N = len(indices)
indices1 = indices[: int(N * 0.7)]
indices2 = indices[int(N * 0.7) :]

# We define a function space for the parameters (at the elements and nodes)

V_dg = dolfin.FunctionSpace(data.mesh, "DG", 0)
V_ode = dolfin.FunctionSpace(data.mesh, "Lagrange", 1)
parameters_ode = np.zeros((len(parameters), V_ode.dim()))
parameters_ode.T[:] = parameters

# And specify a 50% reduction in the `g_K1` conductance

g_K1_index = model["parameter_index"]("g_K1")
# g_K1_index = model["parameter_index"]("scale_drug_IK1")
g_K1_value = parameters[g_K1_index]
g_K1_func = dolfin.Function(V_dg)
g_K1_func.vector()[:] = g_K1_value
g_K1_func.vector()[indices1] = g_K1_value * 0.5
parameters_ode[g_K1_index, :] = (
    dolfin.interpolate(g_K1_func, V_ode).vector().get_local()
)

# and 40% reduction in `g_Na`

g_Na_index = model["parameter_index"]("g_Na")
# g_Na_index = model["parameter_index"]("scale_drug_INa")
g_Na_value = parameters[g_Na_index]
g_Na_func = dolfin.Function(V_dg)
g_Na_func.vector()[:] = g_Na_value
g_Na_func.vector()[indices1] = g_Na_value * 0.6
parameters_ode[g_Na_index, :] = (
    dolfin.interpolate(g_Na_func, V_ode).vector().get_local()
)

# and a 50% reduction in `g_Ca_L`

g_CaL_index = model["parameter_index"]("g_Ca_L")
# g_CaL_index = model["parameter_index"]("scale_drug_ICaL")
g_CaL_value = parameters[g_CaL_index]
g_CaL_func = dolfin.Function(V_dg)
g_CaL_func.vector()[:] = g_CaL_value
g_CaL_func.vector()[indices1] = g_CaL_value * 0.5
parameters_ode[g_CaL_index, :] = (
    dolfin.interpolate(g_CaL_func, V_ode).vector().get_local()
)

# Finally we compute the conductivity tensor with the remaining non conductive tissue

M = define_conductivity_tensor(
    V_dg=V_dg, indices=indices2, f0=data.f0, **conductivities
)

# Now we are ready to set up the models

# +
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
# -

# And solve
# +
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


# ```{figure} ../docs/_static/fibrotic.gif
# ---
# name: reentry
# ---
# ```
