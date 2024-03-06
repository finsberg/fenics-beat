import pytest
import numpy as np
import dolfin

import beat


@pytest.mark.parametrize(
    "M, ac_str, exact, err",
    (
        (
            0.0,
            "cos(2*pi*x[0])*cos(2*pi*x[1])*cos(t)",
            "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)",
            1e-4,
        ),
        (
            1.0,
            "cos(2*pi*x[0])*cos(2*pi*x[1])*(cos(t) +  8*pow(pi, 2)*sin(t))",
            "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)",
            2e-4,
        ),
        (
            2.0,
            "cos(2*pi*x[0])*cos(2*pi*x[1])*(cos(t) +  16*pow(pi, 2)*sin(t))",
            "cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)",
            2e-4,
        ),
    ),
)
def test_monodomain_analtyic(M, ac_str, exact, err):
    N = 15

    time = dolfin.Constant(0.0)
    I_s = dolfin.Expression(ac_str, t=time, degree=5)

    theta = 0.5
    dt = 0.001
    T = 10 * dt

    params = dict(theta=theta, linear_solver_type="direct")  # , default_timestep=dt)

    mesh = dolfin.UnitSquareMesh(N, N)
    model = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)
    res = model.solve((0, T), dt=dt)

    v_exact = dolfin.Expression(exact, t=T, degree=3)

    v_error = dolfin.errornorm(v_exact, res.state, "L2", degree_rise=2)
    assert v_error < err


@pytest.mark.skip_in_parallel
def test_monodomain_spatial_convergence():
    time = dolfin.Constant(0.0)
    ac_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*(cos(t) +  8*pow(pi, 2)*sin(t))"

    I_s = dolfin.Expression(ac_str, t=time, degree=5)

    theta = 0.5
    dt = 0.001
    T = 10 * dt

    params = dict(theta=theta, linear_solver_type="direct")  # , default_timestep=dt)
    v_exact = dolfin.Expression("cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)", t=T, degree=3)

    errors = []
    Ns = [2**level for level in (2, 3, 4, 5)]
    for N in Ns:
        mesh = dolfin.UnitSquareMesh(N, N)
        model = beat.MonodomainModel(
            time=time, mesh=mesh, M=1.0, I_s=I_s, params=params
        )
        res = model.solve((0, T), dt=dt)
        v_error = dolfin.errornorm(v_exact, res.state, "L2", degree_rise=2)
        errors.append(v_error)

    rates = [np.log(e1 / e2) / np.log(2) for e1, e2 in zip(errors[:-1], errors[1:])]
    assert all(rate >= 2.0 for rate in rates)


@pytest.mark.skip_in_parallel
def test_monodomain_temporal_convergence():
    time = dolfin.Constant(0.0)
    ac_str = "cos(2*pi*x[0])*cos(2*pi*x[1])*(cos(t) +  8*pow(pi, 2)*sin(t))"

    I_s = dolfin.Expression(ac_str, t=time, degree=5)

    theta = 0.5
    N = 100
    mesh = dolfin.UnitSquareMesh(N, N)

    T = 1.0

    params = dict(theta=theta, linear_solver_type="direct")  # , default_timestep=dt)
    v_exact = dolfin.Expression("cos(2*pi*x[0])*cos(2*pi*x[1])*sin(t)", t=T, degree=3)

    errors = []
    dts = [1.0 / (2**level) for level in (0, 1, 2, 3)]
    for dt in dts:
        model = beat.MonodomainModel(
            time=time, mesh=mesh, M=1.0, I_s=I_s, params=params
        )
        res = model.solve((0, T), dt=dt)
        v_error = dolfin.errornorm(v_exact, res.state, "L2", degree_rise=2)
        errors.append(v_error)

    rates = [np.log(e1 / e2) / np.log(2) for e1, e2 in zip(errors[:-1], errors[1:])]
    assert all(rate >= 2.0 for rate in rates)
