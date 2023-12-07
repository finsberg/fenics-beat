import pytest
import beat
import numpy as np


@pytest.mark.parametrize(
    "rhs", ("forward_explicit_euler", "forward_generalized_rush_larsen")
)
@pytest.mark.parametrize("model", beat.cellmodels.all_cellmodels())
def test_all_cellmodels(rhs, model):
    states = model.init_state_values()
    old_states = states.copy()

    ode = beat.odesolver.ODESystemSolver(
        fun=getattr(model, rhs),
        states=states,
        parameters=model.init_parameter_values(),
    )

    assert np.allclose(ode.states, old_states)

    ode.step(0.0, 0.01)

    assert not np.allclose(ode.states, old_states)
