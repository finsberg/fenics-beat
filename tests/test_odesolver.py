import numpy as np
import dolfin

import beat
from beat.odesolver import ODESytemSolver, DolfinODESolver


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
            parameters=None,
        )
        j = 0
        t = 0.0
        for _ in range(int((t_bound - t0) / dt)):
            ode.step(t, dt)
            t += dt
            if np.isclose(t, x[j]):
                print(t, j)
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
    t0 = 0.0
    old_states = np.copy(states)

    ode = ODESytemSolver(
        fun=beat.cellmodels.beeler_reuter.forward_generalized_rush_larsen,
        states=states,
        parameters=parameters,
    )
    assert np.allclose(ode.states, old_states)

    ode.step(t0, dt)

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
    t0 = 0.0

    dolfin_ode = DolfinODESolver(
        s,
        fun=beat.cellmodels.beeler_reuter.forward_generalized_rush_larsen,
        init_states=init_states,
        parameters=parameters,
    )
    assert np.allclose(dolfin_ode.s.vector().get_local(), 0.0)
    dolfin_ode.to_dolfin()

    # Just check that values have been updated
    old_state = dolfin_ode.s.vector().get_local().copy()
    assert not np.allclose(old_state, 0.0)

    N = 10
    t = t0
    for _ in range(N):
        dolfin_ode.step(t, dt)
        t += dt

    dolfin_ode.to_dolfin()

    assert not np.allclose(dolfin_ode.s.vector().get_local(), old_state)


def test_assignment_ode():
    model = beat.cellmodels.beeler_reuter
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
    parameters[model.parameter_indices("IstimAmplitude")] = 1.0
    v_index = model.state_indices("V")
    num_states = len(init_states)

    mesh = dolfin.UnitSquareMesh(5, 5)
    V = dolfin.VectorFunctionSpace(mesh, "Lagrange", 1, dim=num_states)
    s = dolfin.Function(V)
    ode = DolfinODESolver(
        s,
        fun=beat.cellmodels.beeler_reuter.forward_generalized_rush_larsen,
        init_states=init_states,
        parameters=parameters,
        v_index=v_index,
    )
    assert np.allclose(ode.v.vector().get_local(), 0)
    assert np.allclose(ode.values[:, 0], init_states)

    ode.v_to_dolfin()
    assert np.allclose(ode.v.vector().get_local(), init_states[v_index])
    # Check that another state is still zero
    assert np.allclose(ode[v_index - 1].vector().get_local(), 0.0)
    ode.to_dolfin()
    assert np.allclose(ode[v_index - 1].vector().get_local(), init_states[v_index - 1])

    # Now update values for v
    ode.values[v_index, :] = 42.0
    assert np.allclose(ode.v.vector().get_local(), init_states[v_index])
    ode.v_to_dolfin()
    assert np.allclose(ode.v.vector().get_local(), 42.0)
    assert np.allclose(ode.s.vector().get_local()[v_index::num_states], 42.0)

    # Now update dolfin function for v
    ode.v.assign(dolfin.Constant(13.0))
    ode.v_from_dolfin()
    assert np.allclose(ode.values[v_index, :], 13.0)
    assert np.allclose(ode.s.vector().get_local()[v_index::num_states], 13.0)
