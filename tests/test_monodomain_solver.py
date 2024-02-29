import dolfin
import numpy as np
import pytest

import beat


def simple_ode_forward_euler(states, t, dt, parameters):
    v, s = states
    values = np.zeros_like(states)
    values[0] = v - s * dt
    values[1] = s + v * dt
    return values


@pytest.mark.parametrize(
    "odespace",
    [
        "CG_1",
        "CG_2",
        "DG_0",
        "DG_1",
        "Quadrature_2",
        "Quadrature_4",
    ],
)
def test_monodomain_splitting_analytic(odespace):
    # Exact solutions
    v_exact_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact_str = "-cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)"

    # Source term
    ac_str = "8*pow(pi, 2)*cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

    family = "Lagrange"
    degree = 1
    N = 50
    mesh = dolfin.UnitSquareMesh(N, N)
    time = dolfin.Constant(0.0)
    I_s = dolfin.Expression(ac_str, t=time, degree=5)
    M = 1.0
    # T = 4.0
    dt = 0.01
    T = dt
    t0 = 0.0

    params = {
        "family": family,
        "degree": degree,
        "linear_solver_type": "direct",
    }

    pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)

    family, degree = odespace.split("_")
    element = dolfin.FiniteElement(
        family, mesh.ufl_cell(), int(degree), quad_scheme="default"
    )
    V_ode = dolfin.FunctionSpace(mesh, element)
    v_ode = dolfin.Function(V_ode)

    s = dolfin.interpolate(s_exact, V_ode)
    s_arr = s.vector().get_local()
    init_states = np.zeros((2, s_arr.size))
    init_states[1, :] = s_arr

    ode = beat.odesolver.DolfinODESolver(
        v_ode=v_ode,
        v_pde=pde.state,
        fun=simple_ode_forward_euler,
        init_states=init_states,
        parameters=None,
        num_states=2,
        v_index=0,
    )
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)
    solver.solve((t0, T), dt=dt)

    v_exact = dolfin.Expression(v_exact_str, t=T, degree=3)
    v_error = dolfin.errornorm(v_exact, pde.state, "L2", degree_rise=2)
    # dolfin.File("v_exact.pvd") << dolfin.interpolate(v_exact, pde.V)
    # dolfin.File("v.pvd") << pde.state

    assert v_error < 0.002


@pytest.mark.parametrize(
    "odespace",
    [
        "CG_1",
        "CG_2",
        "DG_0",
        "DG_1",
        "Quadrature_2",
        "Quadrature_4",
    ],
)
def test_monodomain_splitting_spatial_convergence(odespace):
    print(odespace, "1111111")
    # Exact solutions
    v_exact_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact_str = "-cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)"

    # Source term
    ac_str = "8*pow(pi, 2)*cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"

    s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

    family = "Lagrange"
    degree = 1

    M = 1.0
    dt = 0.01
    T = 1.0
    t0 = 0.0

    params = {
        "family": family,
        "degree": degree,
        "linear_solver_type": "direct",
    }
    # print(f"comm: {mesh.mpi_comm().rank}: 1 - ", N)
    print("a")
    errors = []
    ode_family, ode_degree = odespace.split("_")
    Ns = [2**level for level in (3, 4, 5)]
    for N in Ns:
        mesh = dolfin.UnitSquareMesh(N, N)
        time = dolfin.Constant(0.0)
        I_s = dolfin.Expression(ac_str, t=time, degree=5)

        pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
        s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

        element = dolfin.FiniteElement(
            ode_family, mesh.ufl_cell(), int(ode_degree), quad_scheme="default"
        )
        V_ode = dolfin.FunctionSpace(mesh, element)
        v_ode = dolfin.Function(V_ode)

        s = dolfin.interpolate(s_exact, V_ode)
        s_arr = s.vector().get_local()
        init_states = np.zeros((2, s_arr.size))
        init_states[1, :] = s_arr

        print(f"comm: {mesh.mpi_comm().rank}: 0 - ", init_states.shape)

        ode = beat.odesolver.DolfinODESolver(
            v_ode=v_ode,
            v_pde=pde.state,
            fun=simple_ode_forward_euler,
            init_states=init_states,
            num_states=2,
            parameters=None,
            v_index=0,
        )
        print(f"comm: {mesh.mpi_comm().rank}: 1 - ", N)
        solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)
        solver.solve((t0, T), dt=dt)
        print(f"comm: {mesh.mpi_comm().rank}: 2 - ", N)

        v_exact = dolfin.Expression(v_exact_str, t=T, degree=3)
        v_error = dolfin.errornorm(v_exact, pde.state, "L2", degree_rise=2)
        errors.append(v_error)
        print(f"comm: {mesh.mpi_comm().rank}: 3 - ", N)
        mesh.mpi_comm().barrier()
        print(f"comm: {mesh.mpi_comm().rank}: 4 - ", N)

    rates = [np.log(e1 / e2) / np.log(2) for e1, e2 in zip(errors[:-1], errors[1:])]
    cvg_rate = sum(rates) / len(rates)
    if degree > int(ode_degree):  # DG_0 (ODE) -> CG_1 (PDE)

        assert np.isclose(cvg_rate, int(ode_degree) + 1, rtol=0.1)
    else:
        assert np.isclose(cvg_rate, degree + 1, rtol=0.1)


def test_monodomain_splitting_temporal_convergence():
    # Exact solutions
    v_exact_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact_str = "-cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)"

    # Source term
    ac_str = "8*pow(pi, 2)*cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"

    s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

    family = "Lagrange"
    degree = 1

    M = 1.0
    dt = 0.01
    T = 1.0
    t0 = 0.0

    params = {
        "family": family,
        "degree": degree,
        "linear_solver_type": "direct",
    }

    errors = []
    mesh = dolfin.UnitSquareMesh(150, 150)
    dts = [1.0 / (2**level) for level in (2, 3, 4)]
    for dt in dts:
        time = dolfin.Constant(0.0)
        I_s = dolfin.Expression(ac_str, t=time, degree=5)

        pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
        s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)
        s = dolfin.interpolate(s_exact, pde.V)
        s_arr = s.vector().get_local()
        init_states = np.zeros((2, s_arr.size))
        init_states[1, :] = s_arr

        ode = beat.odesolver.DolfinODESolver(
            v_ode=dolfin.Function(pde.V),
            v_pde=pde.state,
            fun=simple_ode_forward_euler,
            init_states=init_states,
            num_states=2,
            parameters=None,
            v_index=0,
        )
        solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)
        solver.solve((t0, T), dt=dt)

        v_exact = dolfin.Expression(v_exact_str, t=T, degree=3)
        v_error = dolfin.errornorm(v_exact, pde.state, "L2", degree_rise=2)
        errors.append(v_error)

    rates = [np.log(e1 / e2) / np.log(2) for e1, e2 in zip(errors[:-1], errors[1:])]
    cvg_rate = sum(rates) / len(rates)
    # Forward Euler has error of order one in time
    assert np.isclose(cvg_rate, 1, rtol=0.01)


@pytest.mark.parametrize(
    "odespace",
    [
        "CG_1",
        "CG_2",
        "DG_0",
        "DG_1",
        "Quadrature_2",
        "Quadrature_4",
    ],
)
def test_monodomain_splitting_analytic_multiODE(odespace):
    # Exact solutions
    v_exact_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact_str = "-cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)"

    # Source term
    ac_str = "8*pow(pi, 2)*cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

    family = "Lagrange"
    degree = 1
    N = 50
    mesh = dolfin.UnitSquareMesh(N, N)
    time = dolfin.Constant(0.0)
    I_s = dolfin.Expression(ac_str, t=time, degree=5)
    M = 1.0
    # T = 4.0
    dt = 0.01
    T = dt
    t0 = 0.0

    params = {
        "family": family,
        "degree": degree,
        "linear_solver_type": "direct",
    }

    pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)

    family, degree = odespace.split("_")
    element = dolfin.FiniteElement(
        family, mesh.ufl_cell(), int(degree), quad_scheme="default"
    )
    V_ode = dolfin.FunctionSpace(mesh, element)
    v_ode = dolfin.Function(V_ode)

    s = dolfin.interpolate(s_exact, V_ode)
    s_arr = s.vector().get_local()
    init_states = np.zeros((2, s_arr.size))
    init_states[1, :] = s_arr

    ode = beat.odesolver.DolfinODESolver(
        v_ode=v_ode,
        v_pde=pde.state,
        fun=simple_ode_forward_euler,
        init_states=init_states,
        parameters=None,
        num_states=2,
        v_index=0,
    )
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)
    solver.solve((t0, T), dt=dt)

    v_exact = dolfin.Expression(v_exact_str, t=T, degree=3)
    v_error = dolfin.errornorm(v_exact, pde.state, "L2", degree_rise=2)
    # dolfin.File("v_exact.pvd") << dolfin.interpolate(v_exact, pde.V)
    # dolfin.File("v.pvd") << pde.state
    assert v_error < 0.002
