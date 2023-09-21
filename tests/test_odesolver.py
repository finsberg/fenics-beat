import numpy as np
import dolfin

import beat
from beat.odesolver import ODESytemSolver
from beat.dolfin_odesolver import DolfinODESolver


def test_simple_ode_odesystemsolver():
    def simple_ode_forward_euler(states, t, dt, parameters):
        v, s = states
        states[0] = v - s * dt
        states[1] = s + v * dt

    num_points = 1

    t_bound = 1.0
    t0 = 0.0
    x = np.arange(0.1, t_bound + 0.1, 0.1)
    y = np.zeros((len(x), 2))
    sol = np.vstack((np.cos(x), np.sin(x))).T

    errors = []
    for dt in [0.1, 0.01, 0.001, 0.0001]:
        states = np.zeros((2, num_points))
        states.T[:] = [1, 0]
        ode = ODESytemSolver(
            fun=simple_ode_forward_euler,
            states=states,
            dt=dt,
            t_bound=t_bound,
            t0=t0,
            parameters=None,
        )
        j = 0
        for _ in range(int((t_bound - t0) / dt)):
            ode.step()
            if np.isclose(ode.t, x[j]):
                print(ode.t, j)
                y[j, :] = ode.states[:, 0]
                j += 1
        errors.append(np.linalg.norm(sol - y))
    rates = [np.log(e1 / e2) / np.log(10) for e1, e2 in zip(errors[:-1], errors[1:])]
    assert np.allclose(rates, 1, atol=0.01)


def test_beeler_reuter_odesystemsolver():
    model = beat.cellmodels.beeler_reuter
    num_points = 10
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
    parameters[model.parameter_indices("IstimAmplitude")] = 1.0
    num_states = len(init_states)
    states = np.zeros((num_states, num_points))
    states.T[:] = init_states
    dt = 0.1
    t_bound = 1.0
    t0 = 0.0
    old_states = np.copy(states)

    ode = ODESytemSolver(
        fun=beat.cellmodels.beeler_reuter.forward_generalized_rush_larsen,
        states=states,
        dt=dt,
        t_bound=t_bound,
        t0=t0,
        parameters=parameters,
    )
    assert np.isclose(ode.t, t0)
    assert np.allclose(ode.states, old_states)

    ode.step()

    assert np.isclose(ode.t, t0 + dt)
    assert not np.allclose(ode.states, old_states)


def test_beeler_reuter_unit_square():
    model = beat.cellmodels.beeler_reuter
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
    parameters[model.parameter_indices("IstimAmplitude")] = 1.0
    num_states = len(init_states)

    mesh = dolfin.UnitSquareMesh(5, 5)
    V = dolfin.VectorFunctionSpace(mesh, "Lagrange", 1, dim=num_states)
    s = dolfin.Function(V)
    dt = 0.1
    t_bound = 1.0
    t0 = 0.0

    dolfin_ode = DolfinODESolver(
        s,
        fun=beat.cellmodels.beeler_reuter.forward_generalized_rush_larsen,
        init_states=init_states,
        dt=dt,
        t_bound=t_bound,
        t0=t0,
        parameters=parameters,
    )
    assert np.allclose(dolfin_ode.s.vector().get_local(), 0.0)
    dolfin_ode.to_dolfin()

    # Just check that values have been updated
    old_state = dolfin_ode.s.vector().get_local().copy()
    assert not np.allclose(old_state, 0.0)

    N = 10
    for _ in range(N):
        dolfin_ode.step()

    dolfin_ode.to_dolfin()

    assert np.isclose(dolfin_ode.t, N * dt)
    assert not np.allclose(dolfin_ode.s.vector().get_local(), old_state)
