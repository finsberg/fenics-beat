# Gotran generated code for the "tentusscher_panfilov_2006_endo_cell" model

import numpy as np


def init_state_values(**values):
    """
    Initialize state values
    """
    # Init values
    # Xr1=0.00448, Xr2=0.476, Xs=0.0087, m=0.00155, h=0.7573, j=0.7225,
    # d=3.164e-05, f=0.8009, f2=0.9778, fCass=0.9953, s=0.3212,
    # r=2.235e-08, Ca_i=0.00013, R_prime=0.9068, Ca_SR=3.715,
    # Ca_ss=0.00036, Na_i=10.355, V=-86.709, K_i=138.4
    init_values = np.array(
        [
            0.00448,
            0.476,
            0.0087,
            0.00155,
            0.7573,
            0.7225,
            3.164e-05,
            0.8009,
            0.9778,
            0.9953,
            0.3212,
            2.235e-08,
            0.00013,
            0.9068,
            3.715,
            0.00036,
            10.355,
            -86.709,
            138.4,
        ],
        dtype=np.float_,
    )

    # State indices and limit checker
    state_ind = dict(
        [
            ("Xr1", 0),
            ("Xr2", 1),
            ("Xs", 2),
            ("m", 3),
            ("h", 4),
            ("j", 5),
            ("d", 6),
            ("f", 7),
            ("f2", 8),
            ("fCass", 9),
            ("s", 10),
            ("r", 11),
            ("Ca_i", 12),
            ("R_prime", 13),
            ("Ca_SR", 14),
            ("Ca_ss", 15),
            ("Na_i", 16),
            ("V", 17),
            ("K_i", 18),
        ]
    )

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
    # P_kna=0.03, g_K1=5.405, g_Kr=0.153, g_Ks=0.392, g_Na=14.838,
    # g_bna=0.00029, g_CaL=0.0398, g_bca=0.000592, g_to=0.073,
    # K_mNa=40.0, K_mk=1.0, P_NaK=2.724, K_NaCa=1000.0, K_sat=0.1,
    # Km_Ca=1.38, Km_Nai=87.5, alpha=2.5, gamma=0.35, K_pCa=0.0005,
    # g_pCa=0.1238, g_pK=0.0146, Buf_c=0.2, Buf_sr=10.0, Buf_ss=0.4,
    # Ca_o=2.0, EC=1.5, K_buf_c=0.001, K_buf_sr=0.3, K_buf_ss=0.00025,
    # K_up=0.00025, V_leak=0.00036, V_rel=0.102, V_sr=1094.0,
    # V_ss=54.68, V_xfer=0.0038, Vmax_up=0.006375, k1_prime=0.15,
    # k2_prime=0.045, k3=0.06, k4=0.005, max_sr=2.5, min_sr=1,
    # Na_o=140.0, Cm=185.0, F=96.485, R=8.314, T=310.0, V_c=16404.0,
    # stim_amplitude=-52.0, stim_duration=1.0, stim_period=1000.0,
    # stim_start=10.0, K_o=5.4
    init_values = np.array(
        [
            0.03,
            5.405,
            0.153,
            0.392,
            14.838,
            0.00029,
            0.0398,
            0.000592,
            0.073,
            40.0,
            1.0,
            2.724,
            1000.0,
            0.1,
            1.38,
            87.5,
            2.5,
            0.35,
            0.0005,
            0.1238,
            0.0146,
            0.2,
            10.0,
            0.4,
            2.0,
            1.5,
            0.001,
            0.3,
            0.00025,
            0.00025,
            0.00036,
            0.102,
            1094.0,
            54.68,
            0.0038,
            0.006375,
            0.15,
            0.045,
            0.06,
            0.005,
            2.5,
            1,
            140.0,
            185.0,
            96.485,
            8.314,
            310.0,
            16404.0,
            -52.0,
            1.0,
            1000.0,
            10.0,
            5.4,
        ],
        dtype=np.float_,
    )

    # Parameter indices and limit checker
    param_ind = dict(
        [
            ("P_kna", 0),
            ("g_K1", 1),
            ("g_Kr", 2),
            ("g_Ks", 3),
            ("g_Na", 4),
            ("g_bna", 5),
            ("g_CaL", 6),
            ("g_bca", 7),
            ("g_to", 8),
            ("K_mNa", 9),
            ("K_mk", 10),
            ("P_NaK", 11),
            ("K_NaCa", 12),
            ("K_sat", 13),
            ("Km_Ca", 14),
            ("Km_Nai", 15),
            ("alpha", 16),
            ("gamma", 17),
            ("K_pCa", 18),
            ("g_pCa", 19),
            ("g_pK", 20),
            ("Buf_c", 21),
            ("Buf_sr", 22),
            ("Buf_ss", 23),
            ("Ca_o", 24),
            ("EC", 25),
            ("K_buf_c", 26),
            ("K_buf_sr", 27),
            ("K_buf_ss", 28),
            ("K_up", 29),
            ("V_leak", 30),
            ("V_rel", 31),
            ("V_sr", 32),
            ("V_ss", 33),
            ("V_xfer", 34),
            ("Vmax_up", 35),
            ("k1_prime", 36),
            ("k2_prime", 37),
            ("k3", 38),
            ("k4", 39),
            ("max_sr", 40),
            ("min_sr", 41),
            ("Na_o", 42),
            ("Cm", 43),
            ("F", 44),
            ("R", 45),
            ("T", 46),
            ("V_c", 47),
            ("stim_amplitude", 48),
            ("stim_duration", 49),
            ("stim_period", 50),
            ("stim_start", 51),
            ("K_o", 52),
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
    state_inds = dict(
        [
            ("Xr1", 0),
            ("Xr2", 1),
            ("Xs", 2),
            ("m", 3),
            ("h", 4),
            ("j", 5),
            ("d", 6),
            ("f", 7),
            ("f2", 8),
            ("fCass", 9),
            ("s", 10),
            ("r", 11),
            ("Ca_i", 12),
            ("R_prime", 13),
            ("Ca_SR", 14),
            ("Ca_ss", 15),
            ("Na_i", 16),
            ("V", 17),
            ("K_i", 18),
        ]
    )

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
            ("P_kna", 0),
            ("g_K1", 1),
            ("g_Kr", 2),
            ("g_Ks", 3),
            ("g_Na", 4),
            ("g_bna", 5),
            ("g_CaL", 6),
            ("g_bca", 7),
            ("g_to", 8),
            ("K_mNa", 9),
            ("K_mk", 10),
            ("P_NaK", 11),
            ("K_NaCa", 12),
            ("K_sat", 13),
            ("Km_Ca", 14),
            ("Km_Nai", 15),
            ("alpha", 16),
            ("gamma", 17),
            ("K_pCa", 18),
            ("g_pCa", 19),
            ("g_pK", 20),
            ("Buf_c", 21),
            ("Buf_sr", 22),
            ("Buf_ss", 23),
            ("Ca_o", 24),
            ("EC", 25),
            ("K_buf_c", 26),
            ("K_buf_sr", 27),
            ("K_buf_ss", 28),
            ("K_up", 29),
            ("V_leak", 30),
            ("V_rel", 31),
            ("V_sr", 32),
            ("V_ss", 33),
            ("V_xfer", 34),
            ("Vmax_up", 35),
            ("k1_prime", 36),
            ("k2_prime", 37),
            ("k3", 38),
            ("k4", 39),
            ("max_sr", 40),
            ("min_sr", 41),
            ("Na_o", 42),
            ("Cm", 43),
            ("F", 44),
            ("R", 45),
            ("T", 46),
            ("V_c", 47),
            ("stim_amplitude", 48),
            ("stim_duration", 49),
            ("stim_period", 50),
            ("stim_start", 51),
            ("K_o", 52),
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
        [
            ("E_Na", 0),
            ("E_K", 1),
            ("E_Ks", 2),
            ("E_Ca", 3),
            ("alpha_K1", 4),
            ("beta_K1", 5),
            ("xK1_inf", 6),
            ("i_K1", 7),
            ("i_Kr", 8),
            ("xr1_inf", 9),
            ("alpha_xr1", 10),
            ("beta_xr1", 11),
            ("tau_xr1", 12),
            ("xr2_inf", 13),
            ("alpha_xr2", 14),
            ("beta_xr2", 15),
            ("tau_xr2", 16),
            ("i_Ks", 17),
            ("xs_inf", 18),
            ("alpha_xs", 19),
            ("beta_xs", 20),
            ("tau_xs", 21),
            ("i_Na", 22),
            ("m_inf", 23),
            ("alpha_m", 24),
            ("beta_m", 25),
            ("tau_m", 26),
            ("h_inf", 27),
            ("alpha_h", 28),
            ("beta_h", 29),
            ("tau_h", 30),
            ("j_inf", 31),
            ("alpha_j", 32),
            ("beta_j", 33),
            ("tau_j", 34),
            ("i_b_Na", 35),
            ("i_CaL", 36),
            ("d_inf", 37),
            ("alpha_d", 38),
            ("beta_d", 39),
            ("gamma_d", 40),
            ("tau_d", 41),
            ("f_inf", 42),
            ("tau_f", 43),
            ("f2_inf", 44),
            ("tau_f2", 45),
            ("fCass_inf", 46),
            ("tau_fCass", 47),
            ("i_b_Ca", 48),
            ("i_to", 49),
            ("s_inf", 50),
            ("tau_s", 51),
            ("r_inf", 52),
            ("tau_r", 53),
            ("i_NaK", 54),
            ("i_NaCa", 55),
            ("i_p_Ca", 56),
            ("i_p_K", 57),
            ("i_up", 58),
            ("i_leak", 59),
            ("i_xfer", 60),
            ("kcasr", 61),
            ("ddt_Ca_i_total", 62),
            ("f_JCa_i_free", 63),
            ("f_JCa_sr_free", 64),
            ("f_JCa_ss_free", 65),
            ("k1", 66),
            ("k2", 67),
            ("O", 68),
            ("i_rel", 69),
            ("ddt_Ca_sr_total", 70),
            ("ddt_Ca_ss_total", 71),
            ("i_Stim", 72),
            ("dXr1_dt", 73),
            ("dXr2_dt", 74),
            ("dXs_dt", 75),
            ("dm_dt", 76),
            ("dh_dt", 77),
            ("dj_dt", 78),
            ("dd_dt", 79),
            ("df_dt", 80),
            ("df2_dt", 81),
            ("dfCass_dt", 82),
            ("ds_dt", 83),
            ("dr_dt", 84),
            ("dCa_i_dt", 85),
            ("dR_prime_dt", 86),
            ("dCa_SR_dt", 87),
            ("dCa_ss_dt", 88),
            ("dNa_i_dt", 89),
            ("dV_dt", 90),
            ("dK_i_dt", 91),
        ]
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
    Compute the right hand side of the tentusscher_panfilov_2006_endo_cell ODE
    """

    # Assign states
    assert len(states) == 19
    (
        Xr1,
        Xr2,
        Xs,
        m,
        h,
        j,
        d,
        f,
        f2,
        fCass,
        s,
        r,
        Ca_i,
        R_prime,
        Ca_SR,
        Ca_ss,
        Na_i,
        V,
        K_i,
    ) = states

    # Assign parameters
    assert len(parameters) == 53
    (
        P_kna,
        g_K1,
        g_Kr,
        g_Ks,
        g_Na,
        g_bna,
        g_CaL,
        g_bca,
        g_to,
        K_mNa,
        K_mk,
        P_NaK,
        K_NaCa,
        K_sat,
        Km_Ca,
        Km_Nai,
        alpha,
        gamma,
        K_pCa,
        g_pCa,
        g_pK,
        Buf_c,
        Buf_sr,
        Buf_ss,
        Ca_o,
        EC,
        K_buf_c,
        K_buf_sr,
        K_buf_ss,
        K_up,
        V_leak,
        V_rel,
        V_sr,
        V_ss,
        V_xfer,
        Vmax_up,
        k1_prime,
        k2_prime,
        k3,
        k4,
        max_sr,
        min_sr,
        Na_o,
        Cm,
        F,
        R,
        T,
        V_c,
        stim_amplitude,
        stim_duration,
        stim_period,
        stim_start,
        K_o,
    ) = parameters

    # Init return args
    if values is None:
        values = np.zeros((19,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (19,)

    # Expressions for the Reversal potentials component
    E_Na = R * T * np.log(Na_o / Na_i) / F
    E_K = R * T * np.log(K_o / K_i) / F
    E_Ks = R * T * np.log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F
    E_Ca = 0.5 * R * T * np.log(Ca_o / Ca_i) / F

    # Expressions for the Inward rectifier potassium current component
    alpha_K1 = 0.1 / (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
    beta_K1 = (
        0.36787944117144233 * np.exp(0.1 * V - 0.1 * E_K)
        + 3.0606040200802673 * np.exp(0.0002 * V - 0.0002 * E_K)
    ) / (1 + np.exp(0.5 * E_K - 0.5 * V))
    xK1_inf = alpha_K1 / (alpha_K1 + beta_K1)
    i_K1 = 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-E_K + V) * xK1_inf

    # Expressions for the Rapid time dependent potassium current component
    i_Kr = 0.4303314829119352 * g_Kr * np.sqrt(K_o) * (-E_K + V) * Xr1 * Xr2

    # Expressions for the Xr1 gate component
    xr1_inf = 1.0 / (1 + np.exp(-26 / 7 - V / 7))
    alpha_xr1 = 450 / (1 + np.exp(-9 / 2 - V / 10))
    beta_xr1 = 6 / (1 + 13.581324522578193 * np.exp(0.08695652173913043 * V))
    tau_xr1 = alpha_xr1 * beta_xr1
    values[0] = (-Xr1 + xr1_inf) / tau_xr1

    # Expressions for the Xr2 gate component
    xr2_inf = 1.0 / (1 + np.exp(11 / 3 + V / 24))
    alpha_xr2 = 3 / (1 + np.exp(-3 - V / 20))
    beta_xr2 = 1.12 / (1 + np.exp(-3 + V / 20))
    tau_xr2 = alpha_xr2 * beta_xr2
    values[1] = (-Xr2 + xr2_inf) / tau_xr2

    # Expressions for the Slow time dependent potassium current component
    i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V)

    # Expressions for the Xs gate component
    xs_inf = 1.0 / (1 + np.exp(-5 / 14 - V / 14))
    alpha_xs = 1400 / np.sqrt(1 + np.exp(5 / 6 - V / 6))
    beta_xs = 1.0 / (1 + np.exp(-7 / 3 + V / 15))
    tau_xs = 80 + alpha_xs * beta_xs
    values[2] = (-Xs + xs_inf) / tau_xs

    # Expressions for the Fast sodium current component
    i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j

    # Expressions for the m gate component
    m_inf = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
    )
    alpha_m = 1.0 / (1 + np.exp(-12 - V / 5))
    beta_m = 0.1 / (1 + np.exp(7 + V / 5)) + 0.1 / (1 + np.exp(-1 / 4 + V / 200))
    tau_m = alpha_m * beta_m
    values[3] = (-m + m_inf) / tau_m

    # Expressions for the h gate component
    h_inf = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    alpha_h = np.where(V < -40, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V), 0)
    beta_h = np.where(
        V < -40,
        310000 * np.exp(0.3485 * V) + 2.7 * np.exp(0.079 * V),
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V)),
    )
    tau_h = 1.0 / (alpha_h + beta_h)
    values[4] = (-h + h_inf) / tau_h

    # Expressions for the j gate component
    j_inf = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    alpha_j = np.where(
        V < -40,
        (37.78 + V)
        * (-25428 * np.exp(0.2444 * V) - 6.948e-06 * np.exp(-0.04391 * V))
        / (1 + 50262745825.95399 * np.exp(0.311 * V)),
        0,
    )
    beta_j = np.where(
        V < -40,
        0.02424 * np.exp(-0.01052 * V) / (1 + 0.003960868339904256 * np.exp(-0.1378 * V)),
        0.6 * np.exp(0.057 * V) / (1 + 0.040762203978366204 * np.exp(-0.1 * V)),
    )
    tau_j = 1.0 / (alpha_j + beta_j)
    values[5] = (-j + j_inf) / tau_j

    # Expressions for the Sodium background current component
    i_b_Na = g_bna * (-E_Na + V)

    # Expressions for the L_type Ca current component
    i_CaL = (
        4
        * g_CaL
        * (F * F)
        * (-15 + V)
        * (-Ca_o + 0.25 * Ca_ss * np.exp(F * (-30 + 2 * V) / (R * T)))
        * d
        * f
        * f2
        * fCass
        / (R * T * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
    )

    # Expressions for the d gate component
    d_inf = 1.0 / (1 + 0.34415378686541237 * np.exp(-0.13333333333333333 * V))
    alpha_d = 0.25 + 1.4 / (1 + np.exp(-35 / 13 - V / 13))
    beta_d = 1.4 / (1 + np.exp(1 + V / 5))
    gamma_d = 1.0 / (1 + np.exp(5 / 2 - V / 20))
    tau_d = alpha_d * beta_d + gamma_d
    values[6] = (-d + d_inf) / tau_d

    # Expressions for the f gate component
    f_inf = 1.0 / (1 + np.exp(20 / 7 + V / 7))
    tau_f = (
        20
        + 180 / (1 + np.exp(3 + V / 10))
        + 200 / (1 + np.exp(13 / 10 - V / 10))
        + 1102.5 * np.exp(-((27 + V) * (27 + V)) / 225)
    )
    values[7] = (-f + f_inf) / tau_f

    # Expressions for the F2 gate component
    f2_inf = 0.33 + 0.67 / (1 + np.exp(5 + V / 7))
    tau_f2 = (
        31 / (1 + np.exp(5 / 2 - V / 10))
        + 80 / (1 + np.exp(3 + V / 10))
        + 562 * np.exp(-((27 + V) * (27 + V)) / 240)
    )
    values[8] = (-f2 + f2_inf) / tau_f2

    # Expressions for the FCass gate component
    fCass_inf = 0.4 + 0.6 / (1 + 400.0 * (Ca_ss * Ca_ss))
    tau_fCass = 2 + 80 / (1 + 400.0 * (Ca_ss * Ca_ss))
    values[9] = (-fCass + fCass_inf) / tau_fCass

    # Expressions for the Calcium background current component
    i_b_Ca = g_bca * (-E_Ca + V)

    # Expressions for the Transient outward current component
    i_to = g_to * (-E_K + V) * r * s

    # Expressions for the s gate component
    s_inf = 1.0 / (1 + np.exp(28 / 5 + V / 5))
    tau_s = 8 + 1000 * np.exp(-((67 + V) * (67 + V)) / 1000)
    values[10] = (-s + s_inf) / tau_s

    # Expressions for the r gate component
    r_inf = 1.0 / (1 + np.exp(10 / 3 - V / 6))
    tau_r = 0.8 + 9.5 * np.exp(-((40 + V) * (40 + V)) / 1800)
    values[11] = (-r + r_inf) / tau_r

    # Expressions for the Sodium potassium pump current component
    i_NaK = (
        K_o
        * P_NaK
        * Na_i
        / (
            (K_mNa + Na_i)
            * (K_mk + K_o)
            * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
        )
    )

    # Expressions for the Sodium calcium exchanger current component
    i_NaCa = (
        K_NaCa
        * (
            Ca_o * (Na_i * Na_i * Na_i) * np.exp(F * gamma * V / (R * T))
            - alpha * (Na_o * Na_o * Na_o) * Ca_i * np.exp(F * (-1 + gamma) * V / (R * T))
        )
        / (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )

    # Expressions for the Calcium pump current component
    i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i)

    # Expressions for the Potassium pump current component
    i_p_K = g_pK * (-E_K + V) / (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))

    # Expressions for the Calcium dynamics component
    i_up = Vmax_up / (1 + (K_up * K_up) / (Ca_i * Ca_i))
    i_leak = V_leak * (-Ca_i + Ca_SR)
    i_xfer = V_xfer * (-Ca_i + Ca_ss)
    kcasr = max_sr - (max_sr - min_sr) / (1 + (EC * EC) / (Ca_SR * Ca_SR))
    ddt_Ca_i_total = (
        V_sr * (-i_up + i_leak) / V_c
        + Cm * (-i_b_Ca - i_p_Ca + 2 * i_NaCa) / (2 * F * V_c)
        + i_xfer
    )
    f_JCa_i_free = 1.0 / (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
    f_JCa_sr_free = 1.0 / (1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
    f_JCa_ss_free = 1.0 / (1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
    values[12] = ddt_Ca_i_total * f_JCa_i_free
    k1 = k1_prime / kcasr
    k2 = k2_prime * kcasr
    O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1)
    values[13] = k4 * (1 - R_prime) - Ca_ss * R_prime * k2
    i_rel = V_rel * (-Ca_ss + Ca_SR) * O
    ddt_Ca_sr_total = -i_leak - i_rel + i_up
    ddt_Ca_ss_total = V_sr * i_rel / V_ss - V_c * i_xfer / V_ss - Cm * i_CaL / (2 * F * V_ss)
    values[14] = ddt_Ca_sr_total * f_JCa_sr_free
    values[15] = ddt_Ca_ss_total * f_JCa_ss_free

    # Expressions for the Sodium dynamics component
    values[16] = Cm * (-i_Na - i_b_Na - 3 * i_NaCa - 3 * i_NaK) / (F * V_c)

    # Expressions for the Membrane component
    i_Stim = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) >= stim_start,
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
        ),
        stim_amplitude,
        0,
    )
    values[17] = (
        -i_CaL
        - i_K1
        - i_Kr
        - i_Ks
        - i_Na
        - i_NaCa
        - i_NaK
        - i_Stim
        - i_b_Ca
        - i_b_Na
        - i_p_Ca
        - i_p_K
        - i_to
    )

    # Expressions for the Potassium dynamics component
    values[18] = Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + 2 * i_NaK) / (F * V_c)

    # Return results
    return values


def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the tentusscher_panfilov_2006_endo_cell\
        ODE
    """

    # Assign states
    assert len(states) == 19
    (
        Xr1,
        Xr2,
        Xs,
        m,
        h,
        j,
        d,
        f,
        f2,
        fCass,
        s,
        r,
        Ca_i,
        R_prime,
        Ca_SR,
        Ca_ss,
        Na_i,
        V,
        K_i,
    ) = states

    # Assign parameters
    assert len(parameters) == 53
    (
        P_kna,
        g_K1,
        g_Kr,
        g_Ks,
        g_Na,
        g_bna,
        g_CaL,
        g_bca,
        g_to,
        K_mNa,
        K_mk,
        P_NaK,
        K_NaCa,
        K_sat,
        Km_Ca,
        Km_Nai,
        alpha,
        gamma,
        K_pCa,
        g_pCa,
        g_pK,
        Buf_c,
        Buf_sr,
        Buf_ss,
        Ca_o,
        EC,
        K_buf_c,
        K_buf_sr,
        K_buf_ss,
        K_up,
        V_leak,
        V_rel,
        V_sr,
        V_ss,
        V_xfer,
        Vmax_up,
        k1_prime,
        k2_prime,
        k3,
        k4,
        max_sr,
        min_sr,
        Na_o,
        Cm,
        F,
        R,
        T,
        V_c,
        stim_amplitude,
        stim_duration,
        stim_period,
        stim_start,
        K_o,
    ) = parameters

    # Init return args
    if monitored is None:
        monitored = np.zeros((92,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (92,)

    # Expressions for the Reversal potentials component
    monitored[0] = R * T * np.log(Na_o / Na_i) / F
    monitored[1] = R * T * np.log(K_o / K_i) / F
    monitored[2] = R * T * np.log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F
    monitored[3] = 0.5 * R * T * np.log(Ca_o / Ca_i) / F

    # Expressions for the Inward rectifier potassium current component
    monitored[4] = 0.1 / (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * monitored[1]))
    monitored[5] = (
        0.36787944117144233 * np.exp(0.1 * V - 0.1 * monitored[1])
        + 3.0606040200802673 * np.exp(0.0002 * V - 0.0002 * monitored[1])
    ) / (1 + np.exp(0.5 * monitored[1] - 0.5 * V))
    monitored[6] = monitored[4] / (monitored[4] + monitored[5])
    monitored[7] = 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-monitored[1] + V) * monitored[6]

    # Expressions for the Rapid time dependent potassium current component
    monitored[8] = 0.4303314829119352 * g_Kr * np.sqrt(K_o) * (-monitored[1] + V) * Xr1 * Xr2

    # Expressions for the Xr1 gate component
    monitored[9] = 1.0 / (1 + np.exp(-26 / 7 - V / 7))
    monitored[10] = 450 / (1 + np.exp(-9 / 2 - V / 10))
    monitored[11] = 6 / (1 + 13.581324522578193 * np.exp(0.08695652173913043 * V))
    monitored[12] = monitored[10] * monitored[11]
    monitored[73] = (-Xr1 + monitored[9]) / monitored[12]

    # Expressions for the Xr2 gate component
    monitored[13] = 1.0 / (1 + np.exp(11 / 3 + V / 24))
    monitored[14] = 3 / (1 + np.exp(-3 - V / 20))
    monitored[15] = 1.12 / (1 + np.exp(-3 + V / 20))
    monitored[16] = monitored[14] * monitored[15]
    monitored[74] = (-Xr2 + monitored[13]) / monitored[16]

    # Expressions for the Slow time dependent potassium current component
    monitored[17] = g_Ks * (Xs * Xs) * (-monitored[2] + V)

    # Expressions for the Xs gate component
    monitored[18] = 1.0 / (1 + np.exp(-5 / 14 - V / 14))
    monitored[19] = 1400 / np.sqrt(1 + np.exp(5 / 6 - V / 6))
    monitored[20] = 1.0 / (1 + np.exp(-7 / 3 + V / 15))
    monitored[21] = 80 + monitored[19] * monitored[20]
    monitored[75] = (-Xs + monitored[18]) / monitored[21]

    # Expressions for the Fast sodium current component
    monitored[22] = g_Na * (m * m * m) * (-monitored[0] + V) * h * j

    # Expressions for the m gate component
    monitored[23] = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
    )
    monitored[24] = 1.0 / (1 + np.exp(-12 - V / 5))
    monitored[25] = 0.1 / (1 + np.exp(7 + V / 5)) + 0.1 / (1 + np.exp(-1 / 4 + V / 200))
    monitored[26] = monitored[24] * monitored[25]
    monitored[76] = (-m + monitored[23]) / monitored[26]

    # Expressions for the h gate component
    monitored[27] = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    monitored[28] = np.where(V < -40, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V), 0)
    monitored[29] = np.where(
        V < -40,
        310000 * np.exp(0.3485 * V) + 2.7 * np.exp(0.079 * V),
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V)),
    )
    monitored[30] = 1.0 / (monitored[28] + monitored[29])
    monitored[77] = (-h + monitored[27]) / monitored[30]

    # Expressions for the j gate component
    monitored[31] = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    monitored[32] = np.where(
        V < -40,
        (37.78 + V)
        * (-25428 * np.exp(0.2444 * V) - 6.948e-06 * np.exp(-0.04391 * V))
        / (1 + 50262745825.95399 * np.exp(0.311 * V)),
        0,
    )
    monitored[33] = np.where(
        V < -40,
        0.02424 * np.exp(-0.01052 * V) / (1 + 0.003960868339904256 * np.exp(-0.1378 * V)),
        0.6 * np.exp(0.057 * V) / (1 + 0.040762203978366204 * np.exp(-0.1 * V)),
    )
    monitored[34] = 1.0 / (monitored[32] + monitored[33])
    monitored[78] = (-j + monitored[31]) / monitored[34]

    # Expressions for the Sodium background current component
    monitored[35] = g_bna * (-monitored[0] + V)

    # Expressions for the L_type Ca current component
    monitored[36] = (
        4
        * g_CaL
        * (F * F)
        * (-15 + V)
        * (-Ca_o + 0.25 * Ca_ss * np.exp(F * (-30 + 2 * V) / (R * T)))
        * d
        * f
        * f2
        * fCass
        / (R * T * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
    )

    # Expressions for the d gate component
    monitored[37] = 1.0 / (1 + 0.34415378686541237 * np.exp(-0.13333333333333333 * V))
    monitored[38] = 0.25 + 1.4 / (1 + np.exp(-35 / 13 - V / 13))
    monitored[39] = 1.4 / (1 + np.exp(1 + V / 5))
    monitored[40] = 1.0 / (1 + np.exp(5 / 2 - V / 20))
    monitored[41] = monitored[38] * monitored[39] + monitored[40]
    monitored[79] = (-d + monitored[37]) / monitored[41]

    # Expressions for the f gate component
    monitored[42] = 1.0 / (1 + np.exp(20 / 7 + V / 7))
    monitored[43] = (
        20
        + 180 / (1 + np.exp(3 + V / 10))
        + 200 / (1 + np.exp(13 / 10 - V / 10))
        + 1102.5 * np.exp(-((27 + V) * (27 + V)) / 225)
    )
    monitored[80] = (-f + monitored[42]) / monitored[43]

    # Expressions for the F2 gate component
    monitored[44] = 0.33 + 0.67 / (1 + np.exp(5 + V / 7))
    monitored[45] = (
        31 / (1 + np.exp(5 / 2 - V / 10))
        + 80 / (1 + np.exp(3 + V / 10))
        + 562 * np.exp(-((27 + V) * (27 + V)) / 240)
    )
    monitored[81] = (-f2 + monitored[44]) / monitored[45]

    # Expressions for the FCass gate component
    monitored[46] = 0.4 + 0.6 / (1 + 400.0 * (Ca_ss * Ca_ss))
    monitored[47] = 2 + 80 / (1 + 400.0 * (Ca_ss * Ca_ss))
    monitored[82] = (-fCass + monitored[46]) / monitored[47]

    # Expressions for the Calcium background current component
    monitored[48] = g_bca * (-monitored[3] + V)

    # Expressions for the Transient outward current component
    monitored[49] = g_to * (-monitored[1] + V) * r * s

    # Expressions for the s gate component
    monitored[50] = 1.0 / (1 + np.exp(28 / 5 + V / 5))
    monitored[51] = 8 + 1000 * np.exp(-((67 + V) * (67 + V)) / 1000)
    monitored[83] = (-s + monitored[50]) / monitored[51]

    # Expressions for the r gate component
    monitored[52] = 1.0 / (1 + np.exp(10 / 3 - V / 6))
    monitored[53] = 0.8 + 9.5 * np.exp(-((40 + V) * (40 + V)) / 1800)
    monitored[84] = (-r + monitored[52]) / monitored[53]

    # Expressions for the Sodium potassium pump current component
    monitored[54] = (
        K_o
        * P_NaK
        * Na_i
        / (
            (K_mNa + Na_i)
            * (K_mk + K_o)
            * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
        )
    )

    # Expressions for the Sodium calcium exchanger current component
    monitored[55] = (
        K_NaCa
        * (
            Ca_o * (Na_i * Na_i * Na_i) * np.exp(F * gamma * V / (R * T))
            - alpha * (Na_o * Na_o * Na_o) * Ca_i * np.exp(F * (-1 + gamma) * V / (R * T))
        )
        / (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )

    # Expressions for the Calcium pump current component
    monitored[56] = g_pCa * Ca_i / (K_pCa + Ca_i)

    # Expressions for the Potassium pump current component
    monitored[57] = (
        g_pK * (-monitored[1] + V) / (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))
    )

    # Expressions for the Calcium dynamics component
    monitored[58] = Vmax_up / (1 + (K_up * K_up) / (Ca_i * Ca_i))
    monitored[59] = V_leak * (-Ca_i + Ca_SR)
    monitored[60] = V_xfer * (-Ca_i + Ca_ss)
    monitored[61] = max_sr - (max_sr - min_sr) / (1 + (EC * EC) / (Ca_SR * Ca_SR))
    monitored[62] = (
        V_sr * (-monitored[58] + monitored[59]) / V_c
        + Cm * (-monitored[48] - monitored[56] + 2 * monitored[55]) / (2 * F * V_c)
        + monitored[60]
    )
    monitored[63] = 1.0 / (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
    monitored[64] = 1.0 / (1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
    monitored[65] = 1.0 / (1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
    monitored[85] = monitored[62] * monitored[63]
    monitored[66] = k1_prime / monitored[61]
    monitored[67] = k2_prime * monitored[61]
    monitored[68] = (
        (Ca_ss * Ca_ss) * R_prime * monitored[66] / (k3 + (Ca_ss * Ca_ss) * monitored[66])
    )
    monitored[86] = k4 * (1 - R_prime) - Ca_ss * R_prime * monitored[67]
    monitored[69] = V_rel * (-Ca_ss + Ca_SR) * monitored[68]
    monitored[70] = -monitored[59] - monitored[69] + monitored[58]
    monitored[71] = (
        V_sr * monitored[69] / V_ss
        - V_c * monitored[60] / V_ss
        - Cm * monitored[36] / (2 * F * V_ss)
    )
    monitored[87] = monitored[64] * monitored[70]
    monitored[88] = monitored[65] * monitored[71]

    # Expressions for the Sodium dynamics component
    monitored[89] = (
        Cm * (-monitored[22] - monitored[35] - 3 * monitored[54] - 3 * monitored[55]) / (F * V_c)
    )

    # Expressions for the Membrane component
    monitored[72] = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) >= stim_start,
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
        ),
        stim_amplitude,
        0,
    )
    monitored[90] = (
        -monitored[17]
        - monitored[22]
        - monitored[35]
        - monitored[36]
        - monitored[48]
        - monitored[49]
        - monitored[54]
        - monitored[55]
        - monitored[56]
        - monitored[57]
        - monitored[72]
        - monitored[7]
        - monitored[8]
    )

    # Expressions for the Potassium dynamics component
    monitored[91] = (
        Cm
        * (
            -monitored[17]
            - monitored[49]
            - monitored[57]
            - monitored[72]
            - monitored[7]
            - monitored[8]
            + 2 * monitored[54]
        )
        / (F * V_c)
    )

    # Return results
    return monitored


def forward_explicit_euler(states, t, dt, parameters):
    """
    Compute a forward step using the explicit Euler scheme to the\
        tentusscher_panfilov_2006_endo_cell ODE
    """

    # Assign states
    assert len(states) == 19
    (
        Xr1,
        Xr2,
        Xs,
        m,
        h,
        j,
        d,
        f,
        f2,
        fCass,
        s,
        r,
        Ca_i,
        R_prime,
        Ca_SR,
        Ca_ss,
        Na_i,
        V,
        K_i,
    ) = states

    # Assign parameters
    assert len(parameters) == 53
    (
        P_kna,
        g_K1,
        g_Kr,
        g_Ks,
        g_Na,
        g_bna,
        g_CaL,
        g_bca,
        g_to,
        K_mNa,
        K_mk,
        P_NaK,
        K_NaCa,
        K_sat,
        Km_Ca,
        Km_Nai,
        alpha,
        gamma,
        K_pCa,
        g_pCa,
        g_pK,
        Buf_c,
        Buf_sr,
        Buf_ss,
        Ca_o,
        EC,
        K_buf_c,
        K_buf_sr,
        K_buf_ss,
        K_up,
        V_leak,
        V_rel,
        V_sr,
        V_ss,
        V_xfer,
        Vmax_up,
        k1_prime,
        k2_prime,
        k3,
        k4,
        max_sr,
        min_sr,
        Na_o,
        Cm,
        F,
        R,
        T,
        V_c,
        stim_amplitude,
        stim_duration,
        stim_period,
        stim_start,
        K_o,
    ) = parameters

    # Expressions for the Reversal potentials component
    E_Na = R * T * np.log(Na_o / Na_i) / F
    E_K = R * T * np.log(K_o / K_i) / F
    E_Ks = R * T * np.log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F
    E_Ca = 0.5 * R * T * np.log(Ca_o / Ca_i) / F

    # Expressions for the Inward rectifier potassium current component
    alpha_K1 = 0.1 / (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
    beta_K1 = (
        0.36787944117144233 * np.exp(0.1 * V - 0.1 * E_K)
        + 3.0606040200802673 * np.exp(0.0002 * V - 0.0002 * E_K)
    ) / (1 + np.exp(0.5 * E_K - 0.5 * V))
    xK1_inf = alpha_K1 / (alpha_K1 + beta_K1)
    i_K1 = 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-E_K + V) * xK1_inf

    # Expressions for the Rapid time dependent potassium current component
    i_Kr = 0.4303314829119352 * g_Kr * np.sqrt(K_o) * (-E_K + V) * Xr1 * Xr2

    # Expressions for the Xr1 gate component
    xr1_inf = 1.0 / (1 + np.exp(-26 / 7 - V / 7))
    alpha_xr1 = 450 / (1 + np.exp(-9 / 2 - V / 10))
    beta_xr1 = 6 / (1 + 13.581324522578193 * np.exp(0.08695652173913043 * V))
    tau_xr1 = alpha_xr1 * beta_xr1
    dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1
    states[0] = dt * dXr1_dt + Xr1

    # Expressions for the Xr2 gate component
    xr2_inf = 1.0 / (1 + np.exp(11 / 3 + V / 24))
    alpha_xr2 = 3 / (1 + np.exp(-3 - V / 20))
    beta_xr2 = 1.12 / (1 + np.exp(-3 + V / 20))
    tau_xr2 = alpha_xr2 * beta_xr2
    dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2
    states[1] = dt * dXr2_dt + Xr2

    # Expressions for the Slow time dependent potassium current component
    i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V)

    # Expressions for the Xs gate component
    xs_inf = 1.0 / (1 + np.exp(-5 / 14 - V / 14))
    alpha_xs = 1400 / np.sqrt(1 + np.exp(5 / 6 - V / 6))
    beta_xs = 1.0 / (1 + np.exp(-7 / 3 + V / 15))
    tau_xs = 80 + alpha_xs * beta_xs
    dXs_dt = (-Xs + xs_inf) / tau_xs
    states[2] = dt * dXs_dt + Xs

    # Expressions for the Fast sodium current component
    i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j

    # Expressions for the m gate component
    m_inf = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
    )
    alpha_m = 1.0 / (1 + np.exp(-12 - V / 5))
    beta_m = 0.1 / (1 + np.exp(7 + V / 5)) + 0.1 / (1 + np.exp(-1 / 4 + V / 200))
    tau_m = alpha_m * beta_m
    dm_dt = (-m + m_inf) / tau_m
    states[3] = dt * dm_dt + m

    # Expressions for the h gate component
    h_inf = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    alpha_h = np.where(V < -40, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V), 0)
    beta_h = np.where(
        V < -40,
        310000 * np.exp(0.3485 * V) + 2.7 * np.exp(0.079 * V),
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V)),
    )
    tau_h = 1.0 / (alpha_h + beta_h)
    dh_dt = (-h + h_inf) / tau_h
    states[4] = dt * dh_dt + h

    # Expressions for the j gate component
    j_inf = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    alpha_j = np.where(
        V < -40,
        (37.78 + V)
        * (-25428 * np.exp(0.2444 * V) - 6.948e-06 * np.exp(-0.04391 * V))
        / (1 + 50262745825.95399 * np.exp(0.311 * V)),
        0,
    )
    beta_j = np.where(
        V < -40,
        0.02424 * np.exp(-0.01052 * V) / (1 + 0.003960868339904256 * np.exp(-0.1378 * V)),
        0.6 * np.exp(0.057 * V) / (1 + 0.040762203978366204 * np.exp(-0.1 * V)),
    )
    tau_j = 1.0 / (alpha_j + beta_j)
    dj_dt = (-j + j_inf) / tau_j
    states[5] = dt * dj_dt + j

    # Expressions for the Sodium background current component
    i_b_Na = g_bna * (-E_Na + V)

    # Expressions for the L_type Ca current component
    i_CaL = (
        4
        * g_CaL
        * (F * F)
        * (-15 + V)
        * (-Ca_o + 0.25 * Ca_ss * np.exp(F * (-30 + 2 * V) / (R * T)))
        * d
        * f
        * f2
        * fCass
        / (R * T * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
    )

    # Expressions for the d gate component
    d_inf = 1.0 / (1 + 0.34415378686541237 * np.exp(-0.13333333333333333 * V))
    alpha_d = 0.25 + 1.4 / (1 + np.exp(-35 / 13 - V / 13))
    beta_d = 1.4 / (1 + np.exp(1 + V / 5))
    gamma_d = 1.0 / (1 + np.exp(5 / 2 - V / 20))
    tau_d = alpha_d * beta_d + gamma_d
    dd_dt = (-d + d_inf) / tau_d
    states[6] = dt * dd_dt + d

    # Expressions for the f gate component
    f_inf = 1.0 / (1 + np.exp(20 / 7 + V / 7))
    tau_f = (
        20
        + 180 / (1 + np.exp(3 + V / 10))
        + 200 / (1 + np.exp(13 / 10 - V / 10))
        + 1102.5 * np.exp(-((27 + V) * (27 + V)) / 225)
    )
    df_dt = (-f + f_inf) / tau_f
    states[7] = dt * df_dt + f

    # Expressions for the F2 gate component
    f2_inf = 0.33 + 0.67 / (1 + np.exp(5 + V / 7))
    tau_f2 = (
        31 / (1 + np.exp(5 / 2 - V / 10))
        + 80 / (1 + np.exp(3 + V / 10))
        + 562 * np.exp(-((27 + V) * (27 + V)) / 240)
    )
    df2_dt = (-f2 + f2_inf) / tau_f2
    states[8] = dt * df2_dt + f2

    # Expressions for the FCass gate component
    fCass_inf = 0.4 + 0.6 / (1 + 400.0 * (Ca_ss * Ca_ss))
    tau_fCass = 2 + 80 / (1 + 400.0 * (Ca_ss * Ca_ss))
    dfCass_dt = (-fCass + fCass_inf) / tau_fCass
    states[9] = dt * dfCass_dt + fCass

    # Expressions for the Calcium background current component
    i_b_Ca = g_bca * (-E_Ca + V)

    # Expressions for the Transient outward current component
    i_to = g_to * (-E_K + V) * r * s

    # Expressions for the s gate component
    s_inf = 1.0 / (1 + np.exp(28 / 5 + V / 5))
    tau_s = 8 + 1000 * np.exp(-((67 + V) * (67 + V)) / 1000)
    ds_dt = (-s + s_inf) / tau_s
    states[10] = dt * ds_dt + s

    # Expressions for the r gate component
    r_inf = 1.0 / (1 + np.exp(10 / 3 - V / 6))
    tau_r = 0.8 + 9.5 * np.exp(-((40 + V) * (40 + V)) / 1800)
    dr_dt = (-r + r_inf) / tau_r
    states[11] = dt * dr_dt + r

    # Expressions for the Sodium potassium pump current component
    i_NaK = (
        K_o
        * P_NaK
        * Na_i
        / (
            (K_mNa + Na_i)
            * (K_mk + K_o)
            * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
        )
    )

    # Expressions for the Sodium calcium exchanger current component
    i_NaCa = (
        K_NaCa
        * (
            Ca_o * (Na_i * Na_i * Na_i) * np.exp(F * gamma * V / (R * T))
            - alpha * (Na_o * Na_o * Na_o) * Ca_i * np.exp(F * (-1 + gamma) * V / (R * T))
        )
        / (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )

    # Expressions for the Calcium pump current component
    i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i)

    # Expressions for the Potassium pump current component
    i_p_K = g_pK * (-E_K + V) / (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))

    # Expressions for the Calcium dynamics component
    i_up = Vmax_up / (1 + (K_up * K_up) / (Ca_i * Ca_i))
    i_leak = V_leak * (-Ca_i + Ca_SR)
    i_xfer = V_xfer * (-Ca_i + Ca_ss)
    kcasr = max_sr - (max_sr - min_sr) / (1 + (EC * EC) / (Ca_SR * Ca_SR))
    ddt_Ca_i_total = (
        V_sr * (-i_up + i_leak) / V_c
        + Cm * (-i_b_Ca - i_p_Ca + 2 * i_NaCa) / (2 * F * V_c)
        + i_xfer
    )
    f_JCa_i_free = 1.0 / (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
    f_JCa_sr_free = 1.0 / (1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
    f_JCa_ss_free = 1.0 / (1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
    dCa_i_dt = ddt_Ca_i_total * f_JCa_i_free
    states[12] = dt * dCa_i_dt + Ca_i
    k1 = k1_prime / kcasr
    k2 = k2_prime * kcasr
    O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1)
    dR_prime_dt = k4 * (1 - R_prime) - Ca_ss * R_prime * k2
    states[13] = dt * dR_prime_dt + R_prime
    i_rel = V_rel * (-Ca_ss + Ca_SR) * O
    ddt_Ca_sr_total = -i_leak - i_rel + i_up
    ddt_Ca_ss_total = V_sr * i_rel / V_ss - V_c * i_xfer / V_ss - Cm * i_CaL / (2 * F * V_ss)
    dCa_SR_dt = ddt_Ca_sr_total * f_JCa_sr_free
    states[14] = dt * dCa_SR_dt + Ca_SR
    dCa_ss_dt = ddt_Ca_ss_total * f_JCa_ss_free
    states[15] = dt * dCa_ss_dt + Ca_ss

    # Expressions for the Sodium dynamics component
    dNa_i_dt = Cm * (-i_Na - i_b_Na - 3 * i_NaCa - 3 * i_NaK) / (F * V_c)
    states[16] = dt * dNa_i_dt + Na_i

    # Expressions for the Membrane component
    i_Stim = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) >= stim_start,
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
        ),
        stim_amplitude,
        0,
    )
    dV_dt = (
        -i_CaL
        - i_K1
        - i_Kr
        - i_Ks
        - i_Na
        - i_NaCa
        - i_NaK
        - i_Stim
        - i_b_Ca
        - i_b_Na
        - i_p_Ca
        - i_p_K
        - i_to
    )
    states[17] = dt * dV_dt + V

    # Expressions for the Potassium dynamics component
    dK_i_dt = Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + 2 * i_NaK) / (F * V_c)
    states[18] = dt * dK_i_dt + K_i

    # Return results
    return states


def forward_generalized_rush_larsen(states, t, dt, parameters):
    """
    Compute a forward step using the generalised Rush-Larsen (GRL1) scheme to\
        the tentusscher_panfilov_2006_endo_cell ODE
    """

    # Assign states
    assert len(states) == 19
    (
        Xr1,
        Xr2,
        Xs,
        m,
        h,
        j,
        d,
        f,
        f2,
        fCass,
        s,
        r,
        Ca_i,
        R_prime,
        Ca_SR,
        Ca_ss,
        Na_i,
        V,
        K_i,
    ) = states

    # Assign parameters
    assert len(parameters) == 53
    (
        P_kna,
        g_K1,
        g_Kr,
        g_Ks,
        g_Na,
        g_bna,
        g_CaL,
        g_bca,
        g_to,
        K_mNa,
        K_mk,
        P_NaK,
        K_NaCa,
        K_sat,
        Km_Ca,
        Km_Nai,
        alpha,
        gamma,
        K_pCa,
        g_pCa,
        g_pK,
        Buf_c,
        Buf_sr,
        Buf_ss,
        Ca_o,
        EC,
        K_buf_c,
        K_buf_sr,
        K_buf_ss,
        K_up,
        V_leak,
        V_rel,
        V_sr,
        V_ss,
        V_xfer,
        Vmax_up,
        k1_prime,
        k2_prime,
        k3,
        k4,
        max_sr,
        min_sr,
        Na_o,
        Cm,
        F,
        R,
        T,
        V_c,
        stim_amplitude,
        stim_duration,
        stim_period,
        stim_start,
        K_o,
    ) = parameters

    # Expressions for the Reversal potentials component
    E_Na = R * T * np.log(Na_o / Na_i) / F
    E_K = R * T * np.log(K_o / K_i) / F
    E_Ks = R * T * np.log((K_o + Na_o * P_kna) / (P_kna * Na_i + K_i)) / F
    E_Ca = 0.5 * R * T * np.log(Ca_o / Ca_i) / F

    # Expressions for the Inward rectifier potassium current component
    alpha_K1 = 0.1 / (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
    beta_K1 = (
        0.36787944117144233 * np.exp(0.1 * V - 0.1 * E_K)
        + 3.0606040200802673 * np.exp(0.0002 * V - 0.0002 * E_K)
    ) / (1 + np.exp(0.5 * E_K - 0.5 * V))
    xK1_inf = alpha_K1 / (alpha_K1 + beta_K1)
    i_K1 = 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-E_K + V) * xK1_inf

    # Expressions for the Rapid time dependent potassium current component
    i_Kr = 0.4303314829119352 * g_Kr * np.sqrt(K_o) * (-E_K + V) * Xr1 * Xr2

    # Expressions for the Xr1 gate component
    xr1_inf = 1.0 / (1 + np.exp(-26 / 7 - V / 7))
    alpha_xr1 = 450 / (1 + np.exp(-9 / 2 - V / 10))
    beta_xr1 = 6 / (1 + 13.581324522578193 * np.exp(0.08695652173913043 * V))
    tau_xr1 = alpha_xr1 * beta_xr1
    dXr1_dt = (-Xr1 + xr1_inf) / tau_xr1
    dXr1_dt_linearized = -1 / tau_xr1
    states[0] = (-1 + np.exp(dt * dXr1_dt_linearized)) * dXr1_dt / dXr1_dt_linearized + Xr1

    # Expressions for the Xr2 gate component
    xr2_inf = 1.0 / (1 + np.exp(11 / 3 + V / 24))
    alpha_xr2 = 3 / (1 + np.exp(-3 - V / 20))
    beta_xr2 = 1.12 / (1 + np.exp(-3 + V / 20))
    tau_xr2 = alpha_xr2 * beta_xr2
    dXr2_dt = (-Xr2 + xr2_inf) / tau_xr2
    dXr2_dt_linearized = -1 / tau_xr2
    states[1] = (-1 + np.exp(dt * dXr2_dt_linearized)) * dXr2_dt / dXr2_dt_linearized + Xr2

    # Expressions for the Slow time dependent potassium current component
    i_Ks = g_Ks * (Xs * Xs) * (-E_Ks + V)

    # Expressions for the Xs gate component
    xs_inf = 1.0 / (1 + np.exp(-5 / 14 - V / 14))
    alpha_xs = 1400 / np.sqrt(1 + np.exp(5 / 6 - V / 6))
    beta_xs = 1.0 / (1 + np.exp(-7 / 3 + V / 15))
    tau_xs = 80 + alpha_xs * beta_xs
    dXs_dt = (-Xs + xs_inf) / tau_xs
    dXs_dt_linearized = -1 / tau_xs
    states[2] = (-1 + np.exp(dt * dXs_dt_linearized)) * dXs_dt / dXs_dt_linearized + Xs

    # Expressions for the Fast sodium current component
    i_Na = g_Na * (m * m * m) * (-E_Na + V) * h * j

    # Expressions for the m gate component
    m_inf = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * V))
    )
    alpha_m = 1.0 / (1 + np.exp(-12 - V / 5))
    beta_m = 0.1 / (1 + np.exp(7 + V / 5)) + 0.1 / (1 + np.exp(-1 / 4 + V / 200))
    tau_m = alpha_m * beta_m
    dm_dt = (-m + m_inf) / tau_m
    dm_dt_linearized = -1 / tau_m
    states[3] = (-1 + np.exp(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized + m

    # Expressions for the h gate component
    h_inf = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    alpha_h = np.where(V < -40, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * V), 0)
    beta_h = np.where(
        V < -40,
        310000 * np.exp(0.3485 * V) + 2.7 * np.exp(0.079 * V),
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * V)),
    )
    tau_h = 1.0 / (alpha_h + beta_h)
    dh_dt = (-h + h_inf) / tau_h
    dh_dt_linearized = -1 / tau_h
    states[4] = (-1 + np.exp(dt * dh_dt_linearized)) * dh_dt / dh_dt_linearized + h

    # Expressions for the j gate component
    j_inf = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * V))
    )
    alpha_j = np.where(
        V < -40,
        (37.78 + V)
        * (-25428 * np.exp(0.2444 * V) - 6.948e-06 * np.exp(-0.04391 * V))
        / (1 + 50262745825.95399 * np.exp(0.311 * V)),
        0,
    )
    beta_j = np.where(
        V < -40,
        0.02424 * np.exp(-0.01052 * V) / (1 + 0.003960868339904256 * np.exp(-0.1378 * V)),
        0.6 * np.exp(0.057 * V) / (1 + 0.040762203978366204 * np.exp(-0.1 * V)),
    )
    tau_j = 1.0 / (alpha_j + beta_j)
    dj_dt = (-j + j_inf) / tau_j
    dj_dt_linearized = -1 / tau_j
    states[5] = (-1 + np.exp(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized + j

    # Expressions for the Sodium background current component
    i_b_Na = g_bna * (-E_Na + V)

    # Expressions for the L_type Ca current component
    i_CaL = (
        4
        * g_CaL
        * (F * F)
        * (-15 + V)
        * (-Ca_o + 0.25 * Ca_ss * np.exp(F * (-30 + 2 * V) / (R * T)))
        * d
        * f
        * f2
        * fCass
        / (R * T * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
    )

    # Expressions for the d gate component
    d_inf = 1.0 / (1 + 0.34415378686541237 * np.exp(-0.13333333333333333 * V))
    alpha_d = 0.25 + 1.4 / (1 + np.exp(-35 / 13 - V / 13))
    beta_d = 1.4 / (1 + np.exp(1 + V / 5))
    gamma_d = 1.0 / (1 + np.exp(5 / 2 - V / 20))
    tau_d = alpha_d * beta_d + gamma_d
    dd_dt = (-d + d_inf) / tau_d
    dd_dt_linearized = -1 / tau_d
    states[6] = (-1 + np.exp(dt * dd_dt_linearized)) * dd_dt / dd_dt_linearized + d

    # Expressions for the f gate component
    f_inf = 1.0 / (1 + np.exp(20 / 7 + V / 7))
    tau_f = (
        20
        + 180 / (1 + np.exp(3 + V / 10))
        + 200 / (1 + np.exp(13 / 10 - V / 10))
        + 1102.5 * np.exp(-((27 + V) * (27 + V)) / 225)
    )
    df_dt = (-f + f_inf) / tau_f
    df_dt_linearized = -1 / tau_f
    states[7] = (-1 + np.exp(dt * df_dt_linearized)) * df_dt / df_dt_linearized + f

    # Expressions for the F2 gate component
    f2_inf = 0.33 + 0.67 / (1 + np.exp(5 + V / 7))
    tau_f2 = (
        31 / (1 + np.exp(5 / 2 - V / 10))
        + 80 / (1 + np.exp(3 + V / 10))
        + 562 * np.exp(-((27 + V) * (27 + V)) / 240)
    )
    df2_dt = (-f2 + f2_inf) / tau_f2
    df2_dt_linearized = -1 / tau_f2
    states[8] = (-1 + np.exp(dt * df2_dt_linearized)) * df2_dt / df2_dt_linearized + f2

    # Expressions for the FCass gate component
    fCass_inf = 0.4 + 0.6 / (1 + 400.0 * (Ca_ss * Ca_ss))
    tau_fCass = 2 + 80 / (1 + 400.0 * (Ca_ss * Ca_ss))
    dfCass_dt = (-fCass + fCass_inf) / tau_fCass
    dfCass_dt_linearized = -1 / tau_fCass
    states[9] = (-1 + np.exp(dt * dfCass_dt_linearized)) * dfCass_dt / dfCass_dt_linearized + fCass

    # Expressions for the Calcium background current component
    i_b_Ca = g_bca * (-E_Ca + V)

    # Expressions for the Transient outward current component
    i_to = g_to * (-E_K + V) * r * s

    # Expressions for the s gate component
    s_inf = 1.0 / (1 + np.exp(28 / 5 + V / 5))
    tau_s = 8 + 1000 * np.exp(-((67 + V) * (67 + V)) / 1000)
    ds_dt = (-s + s_inf) / tau_s
    ds_dt_linearized = -1 / tau_s
    states[10] = (-1 + np.exp(dt * ds_dt_linearized)) * ds_dt / ds_dt_linearized + s

    # Expressions for the r gate component
    r_inf = 1.0 / (1 + np.exp(10 / 3 - V / 6))
    tau_r = 0.8 + 9.5 * np.exp(-((40 + V) * (40 + V)) / 1800)
    dr_dt = (-r + r_inf) / tau_r
    dr_dt_linearized = -1 / tau_r
    states[11] = (-1 + np.exp(dt * dr_dt_linearized)) * dr_dt / dr_dt_linearized + r

    # Expressions for the Sodium potassium pump current component
    i_NaK = (
        K_o
        * P_NaK
        * Na_i
        / (
            (K_mNa + Na_i)
            * (K_mk + K_o)
            * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
        )
    )

    # Expressions for the Sodium calcium exchanger current component
    i_NaCa = (
        K_NaCa
        * (
            Ca_o * (Na_i * Na_i * Na_i) * np.exp(F * gamma * V / (R * T))
            - alpha * (Na_o * Na_o * Na_o) * Ca_i * np.exp(F * (-1 + gamma) * V / (R * T))
        )
        / (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )

    # Expressions for the Calcium pump current component
    i_p_Ca = g_pCa * Ca_i / (K_pCa + Ca_i)

    # Expressions for the Potassium pump current component
    i_p_K = g_pK * (-E_K + V) / (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))

    # Expressions for the Calcium dynamics component
    i_up = Vmax_up / (1 + (K_up * K_up) / (Ca_i * Ca_i))
    i_leak = V_leak * (-Ca_i + Ca_SR)
    i_xfer = V_xfer * (-Ca_i + Ca_ss)
    kcasr = max_sr - (max_sr - min_sr) / (1 + (EC * EC) / (Ca_SR * Ca_SR))
    ddt_Ca_i_total = (
        V_sr * (-i_up + i_leak) / V_c
        + Cm * (-i_b_Ca - i_p_Ca + 2 * i_NaCa) / (2 * F * V_c)
        + i_xfer
    )
    f_JCa_i_free = 1.0 / (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
    f_JCa_sr_free = 1.0 / (1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
    f_JCa_ss_free = 1.0 / (1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
    dCa_i_dt = ddt_Ca_i_total * f_JCa_i_free
    dE_Ca_dCa_i = -0.5 * R * T / (F * Ca_i)
    dddt_Ca_i_total_di_NaCa = Cm / (F * V_c)
    dddt_Ca_i_total_di_b_Ca = -Cm / (2 * F * V_c)
    dddt_Ca_i_total_di_leak = V_sr / V_c
    dddt_Ca_i_total_di_p_Ca = -Cm / (2 * F * V_c)
    dddt_Ca_i_total_di_up = -V_sr / V_c
    df_JCa_i_free_dCa_i = (
        2
        * Buf_c
        * K_buf_c
        / (
            (
                (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
                * (1 + Buf_c * K_buf_c / ((K_buf_c + Ca_i) * (K_buf_c + Ca_i)))
            )
            * ((K_buf_c + Ca_i) * (K_buf_c + Ca_i) * (K_buf_c + Ca_i))
        )
    )
    di_NaCa_dCa_i = (
        -K_NaCa
        * alpha
        * (Na_o * Na_o * Na_o)
        * np.exp(F * (-1 + gamma) * V / (R * T))
        / (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )
    di_p_Ca_dCa_i = g_pCa / (K_pCa + Ca_i) - g_pCa * Ca_i / ((K_pCa + Ca_i) * (K_pCa + Ca_i))
    di_up_dCa_i = (
        2
        * Vmax_up
        * (K_up * K_up)
        / (
            ((1 + (K_up * K_up) / (Ca_i * Ca_i)) * (1 + (K_up * K_up) / (Ca_i * Ca_i)))
            * (Ca_i * Ca_i * Ca_i)
        )
    )
    dCa_i_dt_linearized = (
        -V_xfer
        + dddt_Ca_i_total_di_NaCa * di_NaCa_dCa_i
        + dddt_Ca_i_total_di_p_Ca * di_p_Ca_dCa_i
        + dddt_Ca_i_total_di_up * di_up_dCa_i
        - V_leak * dddt_Ca_i_total_di_leak
        - g_bca * dE_Ca_dCa_i * dddt_Ca_i_total_di_b_Ca
    ) * f_JCa_i_free + ddt_Ca_i_total * df_JCa_i_free_dCa_i
    states[12] = Ca_i + np.where(
        np.abs(dCa_i_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dCa_i_dt_linearized)) * dCa_i_dt / dCa_i_dt_linearized,
        dt * dCa_i_dt,
    )
    k1 = k1_prime / kcasr
    k2 = k2_prime * kcasr
    O = (Ca_ss * Ca_ss) * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1)
    dR_prime_dt = k4 * (1 - R_prime) - Ca_ss * R_prime * k2
    dR_prime_dt_linearized = -k4 - Ca_ss * k2
    states[13] = (
        np.where(
            np.abs(dR_prime_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dR_prime_dt_linearized)) * dR_prime_dt / dR_prime_dt_linearized,
            dt * dR_prime_dt,
        )
        + R_prime
    )
    i_rel = V_rel * (-Ca_ss + Ca_SR) * O
    ddt_Ca_sr_total = -i_leak - i_rel + i_up
    ddt_Ca_ss_total = V_sr * i_rel / V_ss - V_c * i_xfer / V_ss - Cm * i_CaL / (2 * F * V_ss)
    dCa_SR_dt = ddt_Ca_sr_total * f_JCa_sr_free
    dO_dk1 = (Ca_ss * Ca_ss) * R_prime / (k3 + (Ca_ss * Ca_ss) * k1) - np.power(
        Ca_ss, 4
    ) * R_prime * k1 / ((k3 + (Ca_ss * Ca_ss) * k1) * (k3 + (Ca_ss * Ca_ss) * k1))
    df_JCa_sr_free_dCa_SR = (
        2
        * Buf_sr
        * K_buf_sr
        / (
            (
                (1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
                * (1 + Buf_sr * K_buf_sr / ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR)))
            )
            * ((K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR) * (K_buf_sr + Ca_SR))
        )
    )
    dk1_dkcasr = -k1_prime / (kcasr * kcasr)
    dkcasr_dCa_SR = (
        -2
        * (EC * EC)
        * (max_sr - min_sr)
        / (
            ((1 + (EC * EC) / (Ca_SR * Ca_SR)) * (1 + (EC * EC) / (Ca_SR * Ca_SR)))
            * (Ca_SR * Ca_SR * Ca_SR)
        )
    )
    di_rel_dCa_SR = V_rel * O + V_rel * (-Ca_ss + Ca_SR) * dO_dk1 * dk1_dkcasr * dkcasr_dCa_SR
    di_rel_dO = V_rel * (-Ca_ss + Ca_SR)
    dCa_SR_dt_linearized = (
        -V_leak - di_rel_dCa_SR - dO_dk1 * di_rel_dO * dk1_dkcasr * dkcasr_dCa_SR
    ) * f_JCa_sr_free + ddt_Ca_sr_total * df_JCa_sr_free_dCa_SR
    states[14] = Ca_SR + np.where(
        np.abs(dCa_SR_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dCa_SR_dt_linearized)) * dCa_SR_dt / dCa_SR_dt_linearized,
        dt * dCa_SR_dt,
    )
    dCa_ss_dt = ddt_Ca_ss_total * f_JCa_ss_free
    dO_dCa_ss = -2 * (Ca_ss * Ca_ss * Ca_ss) * (k1 * k1) * R_prime / (
        (k3 + (Ca_ss * Ca_ss) * k1) * (k3 + (Ca_ss * Ca_ss) * k1)
    ) + 2 * Ca_ss * R_prime * k1 / (k3 + (Ca_ss * Ca_ss) * k1)
    dddt_Ca_ss_total_di_CaL = -Cm / (2 * F * V_ss)
    dddt_Ca_ss_total_di_rel = V_sr / V_ss
    dddt_Ca_ss_total_di_xfer = -V_c / V_ss
    df_JCa_ss_free_dCa_ss = (
        2
        * Buf_ss
        * K_buf_ss
        / (
            (
                (1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
                * (1 + Buf_ss * K_buf_ss / ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss)))
            )
            * ((K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss) * (K_buf_ss + Ca_ss))
        )
    )
    di_CaL_dCa_ss = (
        1.0
        * g_CaL
        * (F * F)
        * (-15 + V)
        * d
        * np.exp(F * (-30 + 2 * V) / (R * T))
        * f
        * f2
        * fCass
        / (R * T * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
    )
    di_rel_dCa_ss = -V_rel * O + V_rel * (-Ca_ss + Ca_SR) * dO_dCa_ss
    dCa_ss_dt_linearized = (
        V_xfer * dddt_Ca_ss_total_di_xfer
        + (dO_dCa_ss * di_rel_dO + di_rel_dCa_ss) * dddt_Ca_ss_total_di_rel
        + dddt_Ca_ss_total_di_CaL * di_CaL_dCa_ss
    ) * f_JCa_ss_free + ddt_Ca_ss_total * df_JCa_ss_free_dCa_ss
    states[15] = Ca_ss + np.where(
        np.abs(dCa_ss_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dCa_ss_dt_linearized)) * dCa_ss_dt / dCa_ss_dt_linearized,
        dt * dCa_ss_dt,
    )

    # Expressions for the Sodium dynamics component
    dNa_i_dt = Cm * (-i_Na - i_b_Na - 3 * i_NaCa - 3 * i_NaK) / (F * V_c)
    dE_Na_dNa_i = -R * T / (F * Na_i)
    di_Na_dE_Na = -g_Na * (m * m * m) * h * j
    di_NaCa_dNa_i = (
        3
        * Ca_o
        * K_NaCa
        * (Na_i * Na_i)
        * np.exp(F * gamma * V / (R * T))
        / (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (Ca_o + Km_Ca)
            * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
        )
    )
    di_NaK_dNa_i = K_o * P_NaK / (
        (K_mNa + Na_i)
        * (K_mk + K_o)
        * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
    ) - K_o * P_NaK * Na_i / (
        ((K_mNa + Na_i) * (K_mNa + Na_i))
        * (K_mk + K_o)
        * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
    )
    dNa_i_dt_linearized = (
        Cm
        * (-3 * di_NaCa_dNa_i - 3 * di_NaK_dNa_i + g_bna * dE_Na_dNa_i - dE_Na_dNa_i * di_Na_dE_Na)
        / (F * V_c)
    )
    states[16] = Na_i + np.where(
        np.abs(dNa_i_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dNa_i_dt_linearized)) * dNa_i_dt / dNa_i_dt_linearized,
        dt * dNa_i_dt,
    )

    # Expressions for the Membrane component
    i_Stim = np.where(
        np.logical_and(
            t - stim_period * np.floor(t / stim_period) >= stim_start,
            t - stim_period * np.floor(t / stim_period) <= stim_duration + stim_start,
        ),
        stim_amplitude,
        0,
    )
    dV_dt = (
        -i_CaL
        - i_K1
        - i_Kr
        - i_Ks
        - i_Na
        - i_NaCa
        - i_NaK
        - i_Stim
        - i_b_Ca
        - i_b_Na
        - i_p_Ca
        - i_p_K
        - i_to
    )
    dalpha_K1_dV = (
        -3.686527411996926e-08
        * np.exp(0.06 * V - 0.06 * E_K)
        / (
            (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
            * (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
        )
    )
    dbeta_K1_dV = (
        0.0006121208040160535 * np.exp(0.0002 * V - 0.0002 * E_K)
        + 0.036787944117144235 * np.exp(0.1 * V - 0.1 * E_K)
    ) / (1 + np.exp(0.5 * E_K - 0.5 * V)) + 0.5 * (
        0.36787944117144233 * np.exp(0.1 * V - 0.1 * E_K)
        + 3.0606040200802673 * np.exp(0.0002 * V - 0.0002 * E_K)
    ) * np.exp(0.5 * E_K - 0.5 * V) / (
        (1 + np.exp(0.5 * E_K - 0.5 * V)) * (1 + np.exp(0.5 * E_K - 0.5 * V))
    )
    di_CaL_dV = (
        4
        * g_CaL
        * (F * F)
        * (-Ca_o + 0.25 * Ca_ss * np.exp(F * (-30 + 2 * V) / (R * T)))
        * d
        * f
        * f2
        * fCass
        / (R * T * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
        - 8
        * g_CaL
        * (F * F * F)
        * (-15 + V)
        * (-Ca_o + 0.25 * Ca_ss * np.exp(F * (-30 + 2 * V) / (R * T)))
        * d
        * np.exp(F * (-30 + 2 * V) / (R * T))
        * f
        * f2
        * fCass
        / (
            (R * R)
            * (T * T)
            * (
                (-1 + np.exp(F * (-30 + 2 * V) / (R * T)))
                * (-1 + np.exp(F * (-30 + 2 * V) / (R * T)))
            )
        )
        + 2.0
        * g_CaL
        * (F * F * F)
        * (-15 + V)
        * Ca_ss
        * d
        * np.exp(F * (-30 + 2 * V) / (R * T))
        * f
        * f2
        * fCass
        / ((R * R) * (T * T) * (-1 + np.exp(F * (-30 + 2 * V) / (R * T))))
    )
    dxK1_inf_dalpha_K1 = 1.0 / (alpha_K1 + beta_K1) - alpha_K1 / (
        (alpha_K1 + beta_K1) * (alpha_K1 + beta_K1)
    )
    dxK1_inf_dbeta_K1 = -alpha_K1 / ((alpha_K1 + beta_K1) * (alpha_K1 + beta_K1))
    di_K1_dV = 0.4303314829119352 * g_K1 * np.sqrt(
        K_o
    ) * xK1_inf + 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-E_K + V) * (
        dalpha_K1_dV * dxK1_inf_dalpha_K1 + dbeta_K1_dV * dxK1_inf_dbeta_K1
    )
    di_K1_dxK1_inf = 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-E_K + V)
    di_Kr_dV = 0.4303314829119352 * g_Kr * np.sqrt(K_o) * Xr1 * Xr2
    di_Ks_dV = g_Ks * (Xs * Xs)
    di_Na_dV = g_Na * (m * m * m) * h * j
    di_NaCa_dV = K_NaCa * (
        Ca_o * F * gamma * (Na_i * Na_i * Na_i) * np.exp(F * gamma * V / (R * T)) / (R * T)
        - F
        * alpha
        * (Na_o * Na_o * Na_o)
        * (-1 + gamma)
        * Ca_i
        * np.exp(F * (-1 + gamma) * V / (R * T))
        / (R * T)
    ) / (
        (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
        * (Ca_o + Km_Ca)
        * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
    ) - F * K_NaCa * K_sat * (-1 + gamma) * (
        Ca_o * (Na_i * Na_i * Na_i) * np.exp(F * gamma * V / (R * T))
        - alpha * (Na_o * Na_o * Na_o) * Ca_i * np.exp(F * (-1 + gamma) * V / (R * T))
    ) * np.exp(F * (-1 + gamma) * V / (R * T)) / (
        R
        * T
        * (
            (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
            * (1 + K_sat * np.exp(F * (-1 + gamma) * V / (R * T)))
        )
        * (Ca_o + Km_Ca)
        * ((Km_Nai * Km_Nai * Km_Nai) + (Na_o * Na_o * Na_o))
    )
    di_NaK_dV = (
        K_o
        * P_NaK
        * (
            0.0353 * F * np.exp(-F * V / (R * T)) / (R * T)
            + 0.012450000000000001 * F * np.exp(-0.1 * F * V / (R * T)) / (R * T)
        )
        * Na_i
        / (
            (K_mNa + Na_i)
            * (K_mk + K_o)
            * (
                (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
                * (1 + 0.0353 * np.exp(-F * V / (R * T)) + 0.1245 * np.exp(-0.1 * F * V / (R * T)))
            )
        )
    )
    di_p_K_dV = g_pK / (
        1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V)
    ) + 10.937327047146876 * g_pK * (-E_K + V) * np.exp(-0.16722408026755853 * V) / (
        (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))
        * (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))
    )
    di_to_dV = g_to * r * s
    dV_dt_linearized = (
        -g_bca
        - g_bna
        - di_CaL_dV
        - di_K1_dV
        - di_Kr_dV
        - di_Ks_dV
        - di_NaCa_dV
        - di_NaK_dV
        - di_Na_dV
        - di_p_K_dV
        - di_to_dV
        - (dalpha_K1_dV * dxK1_inf_dalpha_K1 + dbeta_K1_dV * dxK1_inf_dbeta_K1) * di_K1_dxK1_inf
    )
    states[17] = (
        np.where(
            np.abs(dV_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dV_dt_linearized)) * dV_dt / dV_dt_linearized,
            dt * dV_dt,
        )
        + V
    )

    # Expressions for the Potassium dynamics component
    dK_i_dt = Cm * (-i_K1 - i_Kr - i_Ks - i_Stim - i_p_K - i_to + 2 * i_NaK) / (F * V_c)
    dE_K_dK_i = -R * T / (F * K_i)
    dE_Ks_dK_i = -R * T / (F * (P_kna * Na_i + K_i))
    dalpha_K1_dE_K = (
        3.686527411996926e-08
        * np.exp(0.06 * V - 0.06 * E_K)
        / (
            (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
            * (1 + 6.14421235332821e-06 * np.exp(0.06 * V - 0.06 * E_K))
        )
    )
    dbeta_K1_dE_K = (
        -0.0006121208040160535 * np.exp(0.0002 * V - 0.0002 * E_K)
        - 0.036787944117144235 * np.exp(0.1 * V - 0.1 * E_K)
    ) / (1 + np.exp(0.5 * E_K - 0.5 * V)) - 0.5 * (
        0.36787944117144233 * np.exp(0.1 * V - 0.1 * E_K)
        + 3.0606040200802673 * np.exp(0.0002 * V - 0.0002 * E_K)
    ) * np.exp(0.5 * E_K - 0.5 * V) / (
        (1 + np.exp(0.5 * E_K - 0.5 * V)) * (1 + np.exp(0.5 * E_K - 0.5 * V))
    )
    di_K1_dE_K = -0.4303314829119352 * g_K1 * np.sqrt(
        K_o
    ) * xK1_inf + 0.4303314829119352 * g_K1 * np.sqrt(K_o) * (-E_K + V) * (
        dalpha_K1_dE_K * dxK1_inf_dalpha_K1 + dbeta_K1_dE_K * dxK1_inf_dbeta_K1
    )
    di_Kr_dE_K = -0.4303314829119352 * g_Kr * np.sqrt(K_o) * Xr1 * Xr2
    di_Ks_dE_Ks = -g_Ks * (Xs * Xs)
    di_p_K_dE_K = -g_pK / (1 + 65.40521574193832 * np.exp(-0.16722408026755853 * V))
    di_to_dE_K = -g_to * r * s
    dK_i_dt_linearized = (
        Cm
        * (
            -(
                dE_K_dK_i * dalpha_K1_dE_K * dxK1_inf_dalpha_K1
                + dE_K_dK_i * dbeta_K1_dE_K * dxK1_inf_dbeta_K1
            )
            * di_K1_dxK1_inf
            - dE_K_dK_i * di_K1_dE_K
            - dE_K_dK_i * di_Kr_dE_K
            - dE_K_dK_i * di_p_K_dE_K
            - dE_K_dK_i * di_to_dE_K
            - dE_Ks_dK_i * di_Ks_dE_Ks
        )
        / (F * V_c)
    )
    states[18] = K_i + np.where(
        np.abs(dK_i_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dK_i_dt_linearized)) * dK_i_dt / dK_i_dt_linearized,
        dt * dK_i_dt,
    )

    # Return results
    return states
