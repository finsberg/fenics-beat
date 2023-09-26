# Gotran generated code for the "fitzhughnagumo" model

import numpy as np


def init_state_values(**values):
    """
    Initialize state values
    """
    # Init values
    # v=-85.0, s=0.0
    init_values = np.array([-85.0, 0.0], dtype=np.float_)

    # State indices and limit checker
    state_ind = dict([("v", 0), ("s", 1)])

    for state_name, value in values.items():
        if state_name not in state_ind:
            raise ValueError("{0} is not a state.".format(state_name))
        ind = state_ind[state_name]

        # Assign value
        init_values[ind] = value

    return init_values


def init_parameter_values(**values):
    """
    Initialize parameter values
    """
    # Param values
    # a=0.13, b=0.013, c_1=0.26, c_2=0.1, c_3=1.0, stim_amplitude=0,
    # stim_duration=1, stim_period=1000, stim_start=1, v_peak=40.0,
    # v_rest=-85.0
    init_values = np.array(
        [0.13, 0.013, 0.26, 0.1, 1.0, 0, 1, 1000, 1, 40.0, -85.0], dtype=np.float_
    )

    # Parameter indices and limit checker
    param_ind = dict(
        [
            ("a", 0),
            ("b", 1),
            ("c_1", 2),
            ("c_2", 3),
            ("c_3", 4),
            ("stim_amplitude", 5),
            ("stim_duration", 6),
            ("stim_period", 7),
            ("stim_start", 8),
            ("v_peak", 9),
            ("v_rest", 10),
        ]
    )

    for param_name, value in values.items():
        if param_name not in param_ind:
            raise ValueError("{0} is not a parameter.".format(param_name))
        ind = param_ind[param_name]

        # Assign value
        init_values[ind] = value

    return init_values


def state_indices(*states):
    """
    State indices
    """
    state_inds = dict([("v", 0), ("s", 1)])

    indices = []
    for state in states:
        if state not in state_inds:
            raise ValueError("Unknown state: '{0}'".format(state))
        indices.append(state_inds[state])
    if len(indices) > 1:
        return indices
    else:
        return indices[0]


def parameter_indices(*params):
    """
    Parameter indices
    """
    param_inds = dict(
        [
            ("a", 0),
            ("b", 1),
            ("c_1", 2),
            ("c_2", 3),
            ("c_3", 4),
            ("stim_amplitude", 5),
            ("stim_duration", 6),
            ("stim_period", 7),
            ("stim_start", 8),
            ("v_peak", 9),
            ("v_rest", 10),
        ]
    )

    indices = []
    for param in params:
        if param not in param_inds:
            raise ValueError("Unknown param: '{0}'".format(param))
        indices.append(param_inds[param])
    if len(indices) > 1:
        return indices
    else:
        return indices[0]


def monitor_indices(*monitored):
    """
    Monitor indices
    """
    monitor_inds = dict(
        [("v_amp", 0), ("v_th", 1), ("I", 2), ("i_Stim", 3), ("dv_dt", 4), ("ds_dt", 5)]
    )

    indices = []
    for monitor in monitored:
        if monitor not in monitor_inds:
            raise ValueError("Unknown monitored: '{0}'".format(monitor))
        indices.append(monitor_inds[monitor])
    if len(indices) > 1:
        return indices
    else:
        return indices[0]


def rhs(states, t, parameters, values=None):
    """
    Compute the right hand side of the fitzhughnagumo ODE
    """

    # Assign states
    assert len(states) == 2
    v, s = states

    # Assign parameters
    assert len(parameters) == 11
    a = parameters[0]
    b = parameters[1]
    c_1 = parameters[2]
    c_2 = parameters[3]
    c_3 = parameters[4]
    stim_amplitude = parameters[5]
    stim_duration = parameters[6]
    stim_start = parameters[8]
    v_peak = parameters[9]
    v_rest = parameters[10]

    # Init return args
    if values is None:
        values = np.zeros((2,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (2,)

    # Expressions for the fitzhughnagumo component
    v_amp = v_peak - v_rest
    v_th = v_rest + a * v_amp
    I = -c_2 * (-v_rest + v) * s / v_amp + c_1 * (v_peak - v) * (-v_rest + v) * (
        -v_th + v
    ) / (v_amp * v_amp)
    i_Stim = (
        stim_amplitude
        * (1 - 1 / (1 + np.exp(5.0 * t - 5.0 * stim_start)))
        / (1 + np.exp(5.0 * t - 5.0 * stim_duration - 5.0 * stim_start))
    )
    values[0] = I + i_Stim
    values[1] = b * (-v_rest - c_3 * s + v)

    # Return results
    return values


def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the fitzhughnagumo ODE
    """

    # Assign states
    assert len(states) == 2
    v, s = states

    # Assign parameters
    assert len(parameters) == 11
    a = parameters[0]
    b = parameters[1]
    c_1 = parameters[2]
    c_2 = parameters[3]
    c_3 = parameters[4]
    stim_amplitude = parameters[5]
    stim_duration = parameters[6]
    stim_start = parameters[8]
    v_peak = parameters[9]
    v_rest = parameters[10]

    # Init return args
    if monitored is None:
        monitored = np.zeros((6,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (6,)

    # Expressions for the fitzhughnagumo component
    monitored[0] = v_peak - v_rest
    monitored[1] = v_rest + a * monitored[0]
    monitored[2] = -c_2 * (-v_rest + v) * s / monitored[0] + c_1 * (v_peak - v) * (
        -v_rest + v
    ) * (-monitored[1] + v) / (monitored[0] * monitored[0])
    monitored[3] = (
        stim_amplitude
        * (1 - 1 / (1 + np.exp(5.0 * t - 5.0 * stim_start)))
        / (1 + np.exp(5.0 * t - 5.0 * stim_duration - 5.0 * stim_start))
    )
    monitored[4] = monitored[2] + monitored[3]
    monitored[5] = b * (-v_rest - c_3 * s + v)

    # Return results
    return monitored


def forward_explicit_euler(states, t, dt, parameters):
    """
    Compute a forward step using the explicit Euler scheme to the\
        fitzhughnagumo ODE
    """

    # Assign states
    assert len(states) == 2
    v, s = states

    # Assign parameters
    assert len(parameters) == 11
    a = parameters[0]
    b = parameters[1]
    c_1 = parameters[2]
    c_2 = parameters[3]
    c_3 = parameters[4]
    stim_amplitude = parameters[5]
    stim_duration = parameters[6]
    stim_start = parameters[8]
    v_peak = parameters[9]
    v_rest = parameters[10]

    # Expressions for the fitzhughnagumo component
    v_amp = v_peak - v_rest
    v_th = v_rest + a * v_amp
    I = -c_2 * (-v_rest + v) * s / v_amp + c_1 * (v_peak - v) * (-v_rest + v) * (
        -v_th + v
    ) / (v_amp * v_amp)
    i_Stim = (
        stim_amplitude
        * (1 - 1 / (1 + np.exp(5.0 * t - 5.0 * stim_start)))
        / (1 + np.exp(5.0 * t - 5.0 * stim_duration - 5.0 * stim_start))
    )
    dv_dt = I + i_Stim
    states[0] = dt * dv_dt + v
    ds_dt = b * (-v_rest - c_3 * s + v)
    states[1] = dt * ds_dt + s

    # Return results
    return states


def forward_generalized_rush_larsen(states, t, dt, parameters):
    """
    Compute a forward step using the generalised Rush-Larsen (GRL1) scheme to\
        the fitzhughnagumo ODE
    """

    # Assign states
    assert len(states) == 2
    v, s = states

    # Assign parameters
    assert len(parameters) == 11
    a = parameters[0]
    b = parameters[1]
    c_1 = parameters[2]
    c_2 = parameters[3]
    c_3 = parameters[4]
    stim_amplitude = parameters[5]
    stim_duration = parameters[6]
    stim_start = parameters[8]
    v_peak = parameters[9]
    v_rest = parameters[10]

    # Expressions for the fitzhughnagumo component
    v_amp = v_peak - v_rest
    v_th = v_rest + a * v_amp
    I = -c_2 * (-v_rest + v) * s / v_amp + c_1 * (v_peak - v) * (-v_rest + v) * (
        -v_th + v
    ) / (v_amp * v_amp)
    i_Stim = (
        stim_amplitude
        * (1 - 1 / (1 + np.exp(5.0 * t - 5.0 * stim_start)))
        / (1 + np.exp(5.0 * t - 5.0 * stim_duration - 5.0 * stim_start))
    )
    dv_dt = I + i_Stim
    dI_dv = (
        -c_2 * s / v_amp
        + c_1 * (v_peak - v) * (-v_rest + v) / (v_amp * v_amp)
        + c_1 * (v_peak - v) * (-v_th + v) / (v_amp * v_amp)
        - c_1 * (-v_rest + v) * (-v_th + v) / (v_amp * v_amp)
    )
    dv_dt_linearized = dI_dv
    states[0] = (
        np.where(
            np.abs(dv_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dv_dt_linearized)) * dv_dt / dv_dt_linearized,
            dt * dv_dt,
        )
        + v
    )
    ds_dt = b * (-v_rest - c_3 * s + v)
    ds_dt_linearized = -b * c_3
    states[1] = (
        np.where(
            np.abs(ds_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * ds_dt_linearized)) * ds_dt / ds_dt_linearized,
            dt * ds_dt,
        )
        + s
    )

    # Return results
    return states
