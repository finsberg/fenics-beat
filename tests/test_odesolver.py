import numpy as np
import dolfin

import beat
from beat.odesolver import ODESytemSolver
from beat.dolfin_odesolver import DolfinODESolver


def test_beeler_reuter_odesystemsolver():
    model = beat.cellmodels.beeler_reuter
    num_points = 10
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
    # parameters[model.parameter_indices("IstimAmplitude")] = 1.0
    num_states = len(init_states)
    states = np.zeros((num_states, num_points))
    states.T[:] = init_states
    dt = 0.1
    t_bound = 1.0
    t0 = 0.0

    ode = ODESytemSolver(
        fun=beat.cellmodels.beeler_reuter.forward_generalized_rush_larsen,
        states=states,
        dt=dt,
        t_bound=t_bound,
        t0=t0,
        parameters=parameters,
    )
    assert np.isclose(ode.t, t0)
    assert np.allclose(ode.values, states)

    ode.step()

    assert np.isclose(ode.t, t0 + dt)
    assert not np.allclose(ode.values, states)


def test_beeler_reuter_unit_square():
    model = beat.cellmodels.beeler_reuter
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
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
    old_state = dolfin_ode.s.vector().get_local()
    assert not np.allclose(old_state, 0.0)

    dolfin_ode.step()
    dolfin_ode.to_dolfin()

    # parameters[model.parameter_indices("IstimAmplitude")] = 1.0

    # dolfin_ode.values.T[:] = init_states
