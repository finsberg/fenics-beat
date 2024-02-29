# # Niederer benchmark
#
# In this example we will use the same setup as in the Niederer benchmark
# > Niederer SA, Kerfoot E, Benson AP, Bernabeu MO, Bernus O, Bradley C,
#  Cherry EM, Clayton R, Fenton FH, Garny A, Heidenreich E, Land S, Maleckar M,
#  Pathmanathan P, Plank G, Rodríguez JF, Roy I, Sachse FB, Seemann G, Skavhaug O,
#  Smith NP. Verification of cardiac tissue electrophysiology simulators using an
#  N-version benchmark. Philos Trans A Math Phys Eng Sci. 2011 Nov 13;369(1954):4331-51.
#  doi: 10.1098/rsta.2011.0139. PMID: 21969679; PMCID: PMC3263775.
#

import dolfin
import numpy as np
import numpy.typing as npt


import beat

import beat.cellmodels.tentusscher_panfilov_2006.epi as model


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


def setup_geometry(dx):
    Lx = 20.0  # mm
    Ly = 7.0  # mm
    Lz = 3.0  # mm

    mesh = dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
        int(np.rint((Lz / dx))),
    )
    return mesh


def define_stimulus(mesh, chi, C_m, time):
    S1_marker = 1
    L = 1.5
    S1_subdomain = dolfin.CompiledSubDomain(
        "x[0] <= L + DOLFIN_EPS && x[1] <= L + DOLFIN_EPS && x[2] <= L + DOLFIN_EPS",
        L=L,
    )
    S1_markers = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    S1_subdomain.mark(S1_markers, S1_marker)

    # Define stimulation (NB: region of interest carried by the mesh
    # and assumptions in cbcbeat)
    duration = 10.0  # ms
    A = 50000.0  # mu A/cm^3
    cm2mm = 10.0
    factor = 1.0 / (chi * C_m)  # NB: cbcbeat convention
    amplitude = factor * A * (1.0 / cm2mm) ** 3  # mV/ms

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


def define_conductivity_tensor(chi, C_m):
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = 0.17  # mS / mm
    sigma_it = 0.019  # mS / mm
    sigma_el = 0.62  # mS / mm
    sigma_et = 0.24  # mS / mm

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(C_m * chi)
    s_l = sigma_l / (C_m * chi)  # mm^2 / ms
    s_t = sigma_t / (C_m * chi)  # mm^2 / ms

    # Define conductivity tensor
    M = dolfin.as_tensor(((s_l, 0, 0), (0, s_t, 0), (0, 0, s_t)))

    return M


def main():
    fun = model.forward_generalized_rush_larsen
    init_states = setup_initial_conditions()
    parameters = model.init_parameter_values()
    mesh = setup_geometry(dx=0.5)

    # Surface to volume ratio
    chi = 140.0  # mm^{-1}
    # Membrane capacitance
    C_m = 0.01  # mu F / mm^2

    time = dolfin.Constant(0.0)
    I_s = define_stimulus(mesh=mesh, chi=chi, C_m=C_m, time=time)

    M = define_conductivity_tensor(chi, C_m)

    # params = {"linear_solver_type": "direct"}
    params = {"preconditioner": "sor", "use_custom_preconditioner": False}
    ode_space = dolfin.FunctionSpace(mesh, "Lagrange", 1)  # TODO : quadrature space!
    pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
    ode = beat.odesolver.DolfinODESolver(
        v_ode=dolfin.Function(ode_space),
        v_pde=pde.state,
        fun=fun,
        init_states=init_states,
        parameters=parameters,
        num_states=len(init_states),
        v_index=model.state_indices("V"),
    )

    T = 10
    # T = 100  # Change to 100 to reproduce Niederer benchmark
    t = 0.0
    dt = 0.05
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

    xdmf = dolfin.XDMFFile(dolfin.MPI.comm_world, "state.xdmf")
    i = 0
    while t < T + 1e-12:
        if i % 20 == 0:
            v = solver.pde.state.vector().get_local()
            print(f"Solve for {t=:.2f}, {v.max() =}, {v.min() = }")
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


if __name__ == "__main__":
    main()
