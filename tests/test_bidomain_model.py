import dolfin

import beat


def test_bidomain_model():
    dt = 0.001
    T = 10 * dt
    N = 15
    theta = 0.5
    mesh = dolfin.UnitSquareMesh(N, N)
    time = dolfin.Constant(0.0)
    ac_str = "cos(t)*cos(pi*x[0])*cos(pi*x[1]) + pow(pi, 2)*cos(pi*x[0])*cos(pi*x[1])*sin(t)"
    stimulus = dolfin.Expression(ac_str, t=time, degree=5)
    M_i = 1.0
    M_e = 1.0

    # Set-up solver
    params = beat.BidomainModel.default_parameters()
    params["theta"] = theta
    params["linear_solver_type"] = "direct"
    params["use_avg_u_constraint"] = True
    params["enable_adjoint"] = False
    solver = beat.BidomainModel(mesh, time, M_i, M_e, I_s=stimulus, params=params)

    # Define exact solution (Note: v is returned at end of time
    # interval(s), u is computed at somewhere in the time interval
    # depending on theta)
    v_exact = dolfin.Expression("cos(pi*x[0])*cos(pi*x[1])*sin(t)", t=T, degree=3)
    u_exact = dolfin.Expression(
        "-cos(pi*x[0])*cos(pi*x[1])*sin(t)/2.0", t=T - (1.0 - theta) * dt, degree=3
    )

    res = solver.solve((0, T), dt)
    # for timestep, (vs_, vs, vur) in solutions:
    # continue

    v, u, r = res.state.split(deepcopy=True)

    v_error = dolfin.errornorm(v_exact, v, "L2", degree_rise=2)
    u_error = dolfin.errornorm(u_exact, u, "L2", degree_rise=2)
    assert v_error < 2e-5
    assert u_error < 1e-5
