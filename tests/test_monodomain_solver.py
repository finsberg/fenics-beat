import dolfin
import numpy as np
import pytest

import beat


def simple_ode_forward_euler(states, t, dt, parameters):
    v, s = states
    states[0] = v - s * dt
    states[1] = s + v * dt


@pytest.mark.xfail(reason="We have the wrong analytical solution")
def test_monodomain_splitting_analytic():
    # Exact solutions
    v_exact_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact_str = "-cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)"

    # Source term
    ac_str = "8*pow(pi, 2)*cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

    theta = 1.0  # Forward Euler

    family = "Lagrange"
    degree = 1
    N = 50
    mesh = dolfin.UnitSquareMesh(N, N)
    time = dolfin.Constant(0.0)
    I_s = dolfin.Expression(ac_str, t=time, degree=5)
    M = 1.0
    T = 4.0
    dt = 0.01
    t0 = 0.0

    params = {
        "theta": theta,
        "family": family,
        "degree": degree,
        "linear_solver_type": "direct",
    }

    pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)

    s = dolfin.interpolate(s_exact, pde.V)
    s_arr = s.vector().get_local()
    init_states = np.zeros((2, s_arr.size))
    init_states[1, :] = s_arr

    ode = beat.odesolver.DolfinODESolver(
        pde.state,
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
    # breakpoint()
    assert v_error < 1e-4


@pytest.mark.xfail(reason="We have the wrong analytical solution")
def test_monodomain_splitting_spatial_convergence():
    # Exact solutions
    v_exact_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    s_exact_str = "-cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)"

    # Source term
    ac_str = "8*pow(pi, 2)*cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)"
    theta = 1.0  # Forward Euler

    s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)

    family = "Lagrange"
    degree = 1

    M = 1.0
    T = 1.0
    dt = 0.001
    t0 = 0.0

    params = {
        "theta": theta,
        "family": family,
        "degree": degree,
        "linear_solver_type": "direct",
    }

    errors = []
    Ns = [2**level for level in (2, 3, 4, 5)]
    for N in Ns:
        mesh = dolfin.UnitSquareMesh(N, N)
        time = dolfin.Constant(0.0)
        I_s = dolfin.Expression(ac_str, t=time, degree=5)

        pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
        s_exact = dolfin.Expression(s_exact_str, t=0, degree=1)
        s = dolfin.interpolate(s_exact, pde.V)
        s_arr = s.vector().get_local()
        init_states = np.zeros((2, s_arr.size))
        init_states[1, :] = s_arr

        ode = beat.odesolver.DolfinODESolver(
            pde.state,
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
    assert all(rate >= 2.0 for rate in rates)
