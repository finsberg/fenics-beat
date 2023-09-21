import numpy as np

import beat
from beat.cellsolver import ODESytemSolver


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
