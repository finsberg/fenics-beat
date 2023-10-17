from pathlib import Path
import cardiac_geometries

import dolfin
import ufl_legacy as ufl

import beat

import beat.cellmodels.tentusscher_panfilov_2006_epi_cell as model


def get_data(datadir="data_epicardial_stimulation"):
    datadir = Path(datadir)
    msh_file = datadir / "biv_ellipsoid.msh"
    if not msh_file.is_file():
        cardiac_geometries.create_biv_ellipsoid(
            datadir,
            char_length=0.2,  # Reduce this value to get a finer mesh
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


def define_stimulus(mesh, chi, C_m, time, ffun, markers):
    duration = 2.0  # ms
    A = 5  # mu A/cm^3

    factor = 1.0 / (chi * C_m)  # NB: cbcbeat convention
    amplitude = factor * A  # mV/ms

    I_s = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=0.0,
        duration=duration,
        amplitude=amplitude,
        degree=0,
    )

    subdomain_data = dolfin.MeshFunction("size_t", mesh, 2)
    subdomain_data.set_all(0)
    marker = 1
    subdomain_data.array()[ffun.array() == markers["ENDO_LV"][0]] = 1
    subdomain_data.array()[ffun.array() == markers["ENDO_RV"][0]] = 1

    ds = dolfin.Measure("ds", domain=mesh, subdomain_data=subdomain_data)(marker)
    return beat.base_model.Stimulus(dz=ds, expr=I_s)


def define_conductivity_tensor(chi, C_m, f0, s0, n0):
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
    A = dolfin.as_matrix(
        [
            [f0[0], s0[0], n0[0]],
            [f0[1], s0[1], n0[1]],
            [f0[2], s0[2], n0[2]],
        ],
    )

    M_star = ufl.diag(dolfin.as_vector([s_l, s_t, s_t]))
    M = A * M_star * A.T

    return M


def main():
    datadir = Path("data_epicardial_stimulation")
    data = get_data(datadir=datadir)

    fun = model.forward_generalized_rush_larsen
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()

    # Surface to volume ratio
    chi = 140.0  # mm^{-1}
    # Membrane capacitance
    C_m = 0.01  # mu F / mm^2

    time = dolfin.Constant(0.0)
    I_s = define_stimulus(
        mesh=data.mesh,
        chi=chi,
        C_m=C_m,
        time=time,
        ffun=data.ffun,
        markers=data.markers,
    )

    M = define_conductivity_tensor(chi, C_m, f0=data.f0, s0=data.s0, n0=data.n0)

    params = {"preconditioner": "sor", "use_custom_preconditioner": False}
    pde = beat.MonodomainModel(time=time, mesh=data.mesh, M=M, I_s=I_s, params=params)
    ode = beat.odesolver.DolfinODESolver(
        pde.state,
        fun=fun,
        init_states=init_states,
        parameters=parameters,
        num_states=len(init_states),
        v_index=model.state_indices("V"),
    )

    T = 10
    t = 0.0
    dt = 0.05
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

    fname = (datadir / "state.xdmf").as_posix()
    i = 0
    while t < T + 1e-12:
        if i % 20 == 0:
            print(f"Solve for {t=:.2f}, {solver.pde.state.vector().get_local() =}")
            with dolfin.XDMFFile(dolfin.MPI.comm_world, fname) as xdmf:
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
