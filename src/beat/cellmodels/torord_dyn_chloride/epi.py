# Gotran generated code for the "ToRORd_dyn_chloride_epi" model

import numpy as np


def init_state_values(**values):
    """
    Initialize state values
    """
    # Init values
    # CaMKt=0.01273541, d=-2.486527e-36, ff=1.0, fs=0.9510602, fcaf=1.0,
    # fcas=0.9999377, jca=0.9999886, ffp=1.0, fcafp=1.0,
    # nca_ss=0.0003049523, nca_i=0.0005272668, m=0.0005253231,
    # h=0.8645148, j=0.8644571, hp=0.7313656, jp=0.8643527,
    # mL=0.0001117969, hL=0.5916536, hLp=0.3476812, a=0.0008320408,
    # iF=0.9997242, iS=0.9997235, ap=0.0004239121, iFp=0.9997242,
    # iSp=0.9997241, C3=0.0006029079, C2=0.0007393045, C1=0.9984733,
    # O=0.0001787783, I=5.678255e-06, xs1=0.2233584, xs2=0.0001418247,
    # Jrel_np=6.778827e-25, Jrel_p=-1.581941e-23, v=-90.74563,
    # nai=13.40062, nass=13.40094, ki=152.3639, kss=152.3638,
    # cli=34.31721, clss=34.31719, cai=6.621816e-05, cass=5.749921e-05,
    # cansr=1.806794, cajsr=1.805047
    init_values = np.array(
        [
            0.01273541,
            -2.486527e-36,
            1.0,
            0.9510602,
            1.0,
            0.9999377,
            0.9999886,
            1.0,
            1.0,
            0.0003049523,
            0.0005272668,
            0.0005253231,
            0.8645148,
            0.8644571,
            0.7313656,
            0.8643527,
            0.0001117969,
            0.5916536,
            0.3476812,
            0.0008320408,
            0.9997242,
            0.9997235,
            0.0004239121,
            0.9997242,
            0.9997241,
            0.0006029079,
            0.0007393045,
            0.9984733,
            0.0001787783,
            5.678255e-06,
            0.2233584,
            0.0001418247,
            6.778827e-25,
            -1.581941e-23,
            -90.74563,
            13.40062,
            13.40094,
            152.3639,
            152.3638,
            34.31721,
            34.31719,
            6.621816e-05,
            5.749921e-05,
            1.806794,
            1.805047,
        ],
        dtype=np.float_,
    )

    # State indices and limit checker
    state_ind = dict(
        [
            ("CaMKt", 0),
            ("d", 1),
            ("ff", 2),
            ("fs", 3),
            ("fcaf", 4),
            ("fcas", 5),
            ("jca", 6),
            ("ffp", 7),
            ("fcafp", 8),
            ("nca_ss", 9),
            ("nca_i", 10),
            ("m", 11),
            ("h", 12),
            ("j", 13),
            ("hp", 14),
            ("jp", 15),
            ("mL", 16),
            ("hL", 17),
            ("hLp", 18),
            ("a", 19),
            ("iF", 20),
            ("iS", 21),
            ("ap", 22),
            ("iFp", 23),
            ("iSp", 24),
            ("C3", 25),
            ("C2", 26),
            ("C1", 27),
            ("O", 28),
            ("I", 29),
            ("xs1", 30),
            ("xs2", 31),
            ("Jrel_np", 32),
            ("Jrel_p", 33),
            ("v", 34),
            ("nai", 35),
            ("nass", 36),
            ("ki", 37),
            ("kss", 38),
            ("cli", 39),
            ("clss", 40),
            ("cai", 41),
            ("cass", 42),
            ("cansr", 43),
            ("cajsr", 44),
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
    # celltype=1, cao=1.8, clo=150.0, ko=5.0, nao=140.0, F=96485.0, R=8314.0,
    # T=310.0, zca=2, zcl=-1, zk=1, zna=1, L=0.01, rad=0.0011,
    # CaMKo=0.05, KmCaM=0.0015, KmCaMK=0.15, aCaMK=0.05, bCaMK=0.00068,
    # PKNa=0.01833, GNa=11.7802, Gncx_b=0.0034, INaCa_fractionSS=0.35,
    # KmCaAct=0.00015, kasymm=12.5, kcaoff=5000.0, kcaon=1500000.0,
    # kna1=15.0, kna2=5.0, kna3=88.12, qca=0.167, qna=0.5224,
    # wca=60000.0, wna=60000.0, wnaca=5000.0, H=1e-07, Khp=1.698e-07,
    # Kki=0.5, Kko=0.3582, Kmgatp=1.698e-07, Knai0=9.073, Knao0=27.78,
    # Knap=224.0, Kxkur=292.0, MgADP=0.05, MgATP=9.8, Pnak_b=15.4509,
    # delta=-0.155, eP=4.2, k1m=182.4, k1p=949.5, k2m=39.4,
    # k2p=687.2, k3m=79300.0, k3p=1899.0, k4m=40.0, k4p=639.0,
    # PNab=1.9239e-09, GNaL_b=0.0279, thL=200.0, GpCa=0.0005,
    # KmCap=0.0005, tauCa=0.2, tauCl=2.0, tauK=2.0, tauNa=2.0,
    # Aff=0.6, ICaL_fractionSS=0.8, Kmn=0.002, PCa_b=8.3757e-05,
    # dielConstant=74.0, k2n=500.0, offset=0.0, tjca=72.5, vShift=0.0,
    # Jup_b=1.0, A_atp=2.0, K_atp=0.25, K_o_n=5.0, fkatp=0.0,
    # gkatp=4.3195, EKshift=0.0, Gto_b=0.16, GKr_b=0.0321,
    # alpha_1=0.154375, beta_1=0.1911, GKs_b=0.0011, GK1_b=0.6992,
    # GKb_b=0.0189, Fjunc=1, GClCa=0.2843, GClb=0.00198, KdClCa=0.1,
    # PCab=5.9194e-08, Jrel_b=1.5378, bt=4.75, cajsr_half=1.7,
    # i_Stim_Amplitude=-53.0, i_Stim_End=1e+17, i_Stim_Period=1000.0,
    # i_Stim_PulseDuration=1.0, i_Stim_Start=0.0, BSLmax=1.124,
    # BSRmax=0.047, KmBSL=0.0087, KmBSR=0.00087, cmdnmax_b=0.05,
    # csqnmax=10.0, kmcmdn=0.00238, kmcsqn=0.8, kmtrpn=0.0005,
    # trpnmax=0.07
    init_values = np.array(
        [
            1,
            1.8,
            150.0,
            5.0,
            140.0,
            96485.0,
            8314.0,
            310.0,
            2,
            -1,
            1,
            1,
            0.01,
            0.0011,
            0.05,
            0.0015,
            0.15,
            0.05,
            0.00068,
            0.01833,
            11.7802,
            0.0034,
            0.35,
            0.00015,
            12.5,
            5000.0,
            1500000.0,
            15.0,
            5.0,
            88.12,
            0.167,
            0.5224,
            60000.0,
            60000.0,
            5000.0,
            1e-07,
            1.698e-07,
            0.5,
            0.3582,
            1.698e-07,
            9.073,
            27.78,
            224.0,
            292.0,
            0.05,
            9.8,
            15.4509,
            -0.155,
            4.2,
            182.4,
            949.5,
            39.4,
            687.2,
            79300.0,
            1899.0,
            40.0,
            639.0,
            1.9239e-09,
            0.0279,
            200.0,
            0.0005,
            0.0005,
            0.2,
            2.0,
            2.0,
            2.0,
            0.6,
            0.8,
            0.002,
            8.3757e-05,
            74.0,
            500.0,
            0.0,
            72.5,
            0.0,
            1.0,
            2.0,
            0.25,
            5.0,
            0.0,
            4.3195,
            0.0,
            0.16,
            0.0321,
            0.154375,
            0.1911,
            0.0011,
            0.6992,
            0.0189,
            1,
            0.2843,
            0.00198,
            0.1,
            5.9194e-08,
            1.5378,
            4.75,
            1.7,
            -53.0,
            1e17,
            1000.0,
            1.0,
            0.0,
            1.124,
            0.047,
            0.0087,
            0.00087,
            0.05,
            10.0,
            0.00238,
            0.8,
            0.0005,
            0.07,
        ],
        dtype=np.float_,
    )

    # Parameter indices and limit checker
    param_ind = dict(
        [
            ("celltype", 0),
            ("cao", 1),
            ("clo", 2),
            ("ko", 3),
            ("nao", 4),
            ("F", 5),
            ("R", 6),
            ("T", 7),
            ("zca", 8),
            ("zcl", 9),
            ("zk", 10),
            ("zna", 11),
            ("L", 12),
            ("rad", 13),
            ("CaMKo", 14),
            ("KmCaM", 15),
            ("KmCaMK", 16),
            ("aCaMK", 17),
            ("bCaMK", 18),
            ("PKNa", 19),
            ("GNa", 20),
            ("Gncx_b", 21),
            ("INaCa_fractionSS", 22),
            ("KmCaAct", 23),
            ("kasymm", 24),
            ("kcaoff", 25),
            ("kcaon", 26),
            ("kna1", 27),
            ("kna2", 28),
            ("kna3", 29),
            ("qca", 30),
            ("qna", 31),
            ("wca", 32),
            ("wna", 33),
            ("wnaca", 34),
            ("H", 35),
            ("Khp", 36),
            ("Kki", 37),
            ("Kko", 38),
            ("Kmgatp", 39),
            ("Knai0", 40),
            ("Knao0", 41),
            ("Knap", 42),
            ("Kxkur", 43),
            ("MgADP", 44),
            ("MgATP", 45),
            ("Pnak_b", 46),
            ("delta", 47),
            ("eP", 48),
            ("k1m", 49),
            ("k1p", 50),
            ("k2m", 51),
            ("k2p", 52),
            ("k3m", 53),
            ("k3p", 54),
            ("k4m", 55),
            ("k4p", 56),
            ("PNab", 57),
            ("GNaL_b", 58),
            ("thL", 59),
            ("GpCa", 60),
            ("KmCap", 61),
            ("tauCa", 62),
            ("tauCl", 63),
            ("tauK", 64),
            ("tauNa", 65),
            ("Aff", 66),
            ("ICaL_fractionSS", 67),
            ("Kmn", 68),
            ("PCa_b", 69),
            ("dielConstant", 70),
            ("k2n", 71),
            ("offset", 72),
            ("tjca", 73),
            ("vShift", 74),
            ("Jup_b", 75),
            ("A_atp", 76),
            ("K_atp", 77),
            ("K_o_n", 78),
            ("fkatp", 79),
            ("gkatp", 80),
            ("EKshift", 81),
            ("Gto_b", 82),
            ("GKr_b", 83),
            ("alpha_1", 84),
            ("beta_1", 85),
            ("GKs_b", 86),
            ("GK1_b", 87),
            ("GKb_b", 88),
            ("Fjunc", 89),
            ("GClCa", 90),
            ("GClb", 91),
            ("KdClCa", 92),
            ("PCab", 93),
            ("Jrel_b", 94),
            ("bt", 95),
            ("cajsr_half", 96),
            ("i_Stim_Amplitude", 97),
            ("i_Stim_End", 98),
            ("i_Stim_Period", 99),
            ("i_Stim_PulseDuration", 100),
            ("i_Stim_Start", 101),
            ("BSLmax", 102),
            ("BSRmax", 103),
            ("KmBSL", 104),
            ("KmBSR", 105),
            ("cmdnmax_b", 106),
            ("csqnmax", 107),
            ("kmcmdn", 108),
            ("kmcsqn", 109),
            ("kmtrpn", 110),
            ("trpnmax", 111),
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
            ("CaMKt", 0),
            ("d", 1),
            ("ff", 2),
            ("fs", 3),
            ("fcaf", 4),
            ("fcas", 5),
            ("jca", 6),
            ("ffp", 7),
            ("fcafp", 8),
            ("nca_ss", 9),
            ("nca_i", 10),
            ("m", 11),
            ("h", 12),
            ("j", 13),
            ("hp", 14),
            ("jp", 15),
            ("mL", 16),
            ("hL", 17),
            ("hLp", 18),
            ("a", 19),
            ("iF", 20),
            ("iS", 21),
            ("ap", 22),
            ("iFp", 23),
            ("iSp", 24),
            ("C3", 25),
            ("C2", 26),
            ("C1", 27),
            ("O", 28),
            ("I", 29),
            ("xs1", 30),
            ("xs2", 31),
            ("Jrel_np", 32),
            ("Jrel_p", 33),
            ("v", 34),
            ("nai", 35),
            ("nass", 36),
            ("ki", 37),
            ("kss", 38),
            ("cli", 39),
            ("clss", 40),
            ("cai", 41),
            ("cass", 42),
            ("cansr", 43),
            ("cajsr", 44),
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
            ("celltype", 0),
            ("cao", 1),
            ("clo", 2),
            ("ko", 3),
            ("nao", 4),
            ("F", 5),
            ("R", 6),
            ("T", 7),
            ("zca", 8),
            ("zcl", 9),
            ("zk", 10),
            ("zna", 11),
            ("L", 12),
            ("rad", 13),
            ("CaMKo", 14),
            ("KmCaM", 15),
            ("KmCaMK", 16),
            ("aCaMK", 17),
            ("bCaMK", 18),
            ("PKNa", 19),
            ("GNa", 20),
            ("Gncx_b", 21),
            ("INaCa_fractionSS", 22),
            ("KmCaAct", 23),
            ("kasymm", 24),
            ("kcaoff", 25),
            ("kcaon", 26),
            ("kna1", 27),
            ("kna2", 28),
            ("kna3", 29),
            ("qca", 30),
            ("qna", 31),
            ("wca", 32),
            ("wna", 33),
            ("wnaca", 34),
            ("H", 35),
            ("Khp", 36),
            ("Kki", 37),
            ("Kko", 38),
            ("Kmgatp", 39),
            ("Knai0", 40),
            ("Knao0", 41),
            ("Knap", 42),
            ("Kxkur", 43),
            ("MgADP", 44),
            ("MgATP", 45),
            ("Pnak_b", 46),
            ("delta", 47),
            ("eP", 48),
            ("k1m", 49),
            ("k1p", 50),
            ("k2m", 51),
            ("k2p", 52),
            ("k3m", 53),
            ("k3p", 54),
            ("k4m", 55),
            ("k4p", 56),
            ("PNab", 57),
            ("GNaL_b", 58),
            ("thL", 59),
            ("GpCa", 60),
            ("KmCap", 61),
            ("tauCa", 62),
            ("tauCl", 63),
            ("tauK", 64),
            ("tauNa", 65),
            ("Aff", 66),
            ("ICaL_fractionSS", 67),
            ("Kmn", 68),
            ("PCa_b", 69),
            ("dielConstant", 70),
            ("k2n", 71),
            ("offset", 72),
            ("tjca", 73),
            ("vShift", 74),
            ("Jup_b", 75),
            ("A_atp", 76),
            ("K_atp", 77),
            ("K_o_n", 78),
            ("fkatp", 79),
            ("gkatp", 80),
            ("EKshift", 81),
            ("Gto_b", 82),
            ("GKr_b", 83),
            ("alpha_1", 84),
            ("beta_1", 85),
            ("GKs_b", 86),
            ("GK1_b", 87),
            ("GKb_b", 88),
            ("Fjunc", 89),
            ("GClCa", 90),
            ("GClb", 91),
            ("KdClCa", 92),
            ("PCab", 93),
            ("Jrel_b", 94),
            ("bt", 95),
            ("cajsr_half", 96),
            ("i_Stim_Amplitude", 97),
            ("i_Stim_End", 98),
            ("i_Stim_Period", 99),
            ("i_Stim_PulseDuration", 100),
            ("i_Stim_Start", 101),
            ("BSLmax", 102),
            ("BSRmax", 103),
            ("KmBSL", 104),
            ("KmBSR", 105),
            ("cmdnmax_b", 106),
            ("csqnmax", 107),
            ("kmcmdn", 108),
            ("kmcsqn", 109),
            ("kmtrpn", 110),
            ("trpnmax", 111),
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
            ("vcell", 0),
            ("Ageo", 1),
            ("Acap", 2),
            ("vmyo", 3),
            ("vnsr", 4),
            ("vjsr", 5),
            ("vss", 6),
            ("CaMKb", 7),
            ("CaMKa", 8),
            ("ENa", 9),
            ("EK", 10),
            ("EKs", 11),
            ("ECl", 12),
            ("EClss", 13),
            ("mss", 14),
            ("tm", 15),
            ("hss", 16),
            ("ah", 17),
            ("bh", 18),
            ("th", 19),
            ("aj", 20),
            ("bj", 21),
            ("jss", 22),
            ("tj", 23),
            ("hssp", 24),
            ("tjp", 25),
            ("fINap", 26),
            ("INa", 27),
            ("hca", 28),
            ("hna", 29),
            ("h1_i", 30),
            ("h2_i", 31),
            ("h3_i", 32),
            ("h4_i", 33),
            ("h5_i", 34),
            ("h6_i", 35),
            ("h7_i", 36),
            ("h8_i", 37),
            ("h9_i", 38),
            ("h10_i", 39),
            ("h11_i", 40),
            ("h12_i", 41),
            ("k1_i", 42),
            ("k2_i", 43),
            ("k3p_i", 44),
            ("k3pp_i", 45),
            ("k3_i", 46),
            ("k4p_i", 47),
            ("k4pp_i", 48),
            ("k4_i", 49),
            ("k5_i", 50),
            ("k6_i", 51),
            ("k7_i", 52),
            ("k8_i", 53),
            ("x1_i", 54),
            ("x2_i", 55),
            ("x3_i", 56),
            ("x4_i", 57),
            ("E1_i", 58),
            ("E2_i", 59),
            ("E3_i", 60),
            ("E4_i", 61),
            ("allo_i", 62),
            ("JncxNa_i", 63),
            ("JncxCa_i", 64),
            ("Gncx", 65),
            ("INaCa_i", 66),
            ("h1_ss", 67),
            ("h2_ss", 68),
            ("h3_ss", 69),
            ("h4_ss", 70),
            ("h5_ss", 71),
            ("h6_ss", 72),
            ("h7_ss", 73),
            ("h8_ss", 74),
            ("h9_ss", 75),
            ("h10_ss", 76),
            ("h11_ss", 77),
            ("h12_ss", 78),
            ("k1_ss", 79),
            ("k2_ss", 80),
            ("k3p_ss", 81),
            ("k3pp_ss", 82),
            ("k3_ss", 83),
            ("k4p_ss", 84),
            ("k4pp_ss", 85),
            ("k4_ss", 86),
            ("k5_ss", 87),
            ("k6_ss", 88),
            ("k7_ss", 89),
            ("k8_ss", 90),
            ("x1_ss", 91),
            ("x2_ss", 92),
            ("x3_ss", 93),
            ("x4_ss", 94),
            ("E1_ss", 95),
            ("E2_ss", 96),
            ("E3_ss", 97),
            ("E4_ss", 98),
            ("allo_ss", 99),
            ("JncxNa_ss", 100),
            ("JncxCa_ss", 101),
            ("INaCa_ss", 102),
            ("Knai", 103),
            ("Knao", 104),
            ("P", 105),
            ("a1", 106),
            ("b1", 107),
            ("a2", 108),
            ("b2", 109),
            ("a3", 110),
            ("b3", 111),
            ("a4", 112),
            ("b4", 113),
            ("x1", 114),
            ("x2", 115),
            ("x3", 116),
            ("x4", 117),
            ("E1", 118),
            ("E2", 119),
            ("E3", 120),
            ("E4", 121),
            ("JnakNa", 122),
            ("JnakK", 123),
            ("Pnak", 124),
            ("INaK", 125),
            ("INab", 126),
            ("mLss", 127),
            ("tmL", 128),
            ("hLss", 129),
            ("hLssp", 130),
            ("thLp", 131),
            ("GNaL", 132),
            ("fINaLp", 133),
            ("INaL", 134),
            ("IpCa", 135),
            ("JdiffNa", 136),
            ("JdiffK", 137),
            ("Jdiff", 138),
            ("JdiffCl", 139),
            ("dss", 140),
            ("td", 141),
            ("fss", 142),
            ("tff", 143),
            ("tfs", 144),
            ("Afs", 145),
            ("f", 146),
            ("fcass", 147),
            ("tfcaf", 148),
            ("tfcas", 149),
            ("Afcaf", 150),
            ("Afcas", 151),
            ("fca", 152),
            ("jcass", 153),
            ("tffp", 154),
            ("fp", 155),
            ("tfcafp", 156),
            ("fcap", 157),
            ("km2n", 158),
            ("anca_ss", 159),
            ("Io", 160),
            ("Iss", 161),
            ("constA", 162),
            ("gamma_cass", 163),
            ("gamma_cao", 164),
            ("gamma_nass", 165),
            ("gamma_nao", 166),
            ("gamma_kss", 167),
            ("gamma_ko", 168),
            ("PhiCaL_ss", 169),
            ("PhiCaNa_ss", 170),
            ("PhiCaK_ss", 171),
            ("PCa", 172),
            ("PCap", 173),
            ("PCaNa", 174),
            ("PCaK", 175),
            ("PCaNap", 176),
            ("PCaKp", 177),
            ("fICaLp", 178),
            ("ICaL_ss", 179),
            ("ICaNa_ss", 180),
            ("ICaK_ss", 181),
            ("anca_i", 182),
            ("Ii", 183),
            ("gamma_cai", 184),
            ("gamma_nai", 185),
            ("gamma_ki", 186),
            ("PhiCaL_i", 187),
            ("PhiCaNa_i", 188),
            ("PhiCaK_i", 189),
            ("ICaL_i", 190),
            ("ICaNa_i", 191),
            ("ICaK_i", 192),
            ("ICaL", 193),
            ("ICaNa", 194),
            ("ICaK", 195),
            ("upScale", 196),
            ("Jupnp", 197),
            ("Jupp", 198),
            ("fJupp", 199),
            ("Jleak", 200),
            ("Jup", 201),
            ("akik", 202),
            ("bkik", 203),
            ("I_katp", 204),
            ("ass", 205),
            ("ta", 206),
            ("iss", 207),
            ("delta_epi", 208),
            ("tiF_b", 209),
            ("tiS_b", 210),
            ("tiF", 211),
            ("tiS", 212),
            ("AiF", 213),
            ("AiS", 214),
            ("i", 215),
            ("assp", 216),
            ("dti_develop", 217),
            ("dti_recover", 218),
            ("tiFp", 219),
            ("tiSp", 220),
            ("ip", 221),
            ("Gto", 222),
            ("fItop", 223),
            ("Ito", 224),
            ("alpha", 225),
            ("beta", 226),
            ("alpha_2", 227),
            ("beta_2", 228),
            ("alpha_i", 229),
            ("beta_i", 230),
            ("alpha_C2ToI", 231),
            ("beta_ItoC2", 232),
            ("GKr", 233),
            ("IKr", 234),
            ("xs1ss", 235),
            ("txs1", 236),
            ("xs2ss", 237),
            ("txs2", 238),
            ("KsCa", 239),
            ("GKs", 240),
            ("IKs", 241),
            ("aK1", 242),
            ("bK1", 243),
            ("K1ss", 244),
            ("GK1", 245),
            ("IK1", 246),
            ("xkb", 247),
            ("GKb", 248),
            ("IKb", 249),
            ("IClCa_junc", 250),
            ("IClCa_sl", 251),
            ("IClCa", 252),
            ("IClb", 253),
            ("ICab", 254),
            ("a_rel", 255),
            ("Jrel_inf_b", 256),
            ("Jrel_inf", 257),
            ("tau_rel_b", 258),
            ("tau_rel", 259),
            ("btp", 260),
            ("a_relp", 261),
            ("Jrel_infp_b", 262),
            ("Jrel_infp", 263),
            ("tau_relp_b", 264),
            ("tau_relp", 265),
            ("fJrelp", 266),
            ("Jrel", 267),
            ("Istim", 268),
            ("cmdnmax", 269),
            ("Bcai", 270),
            ("Bcass", 271),
            ("Bcajsr", 272),
            ("vfrt", 273),
            ("vffrt", 274),
            ("Jtr", 275),
            ("dCaMKt_dt", 276),
            ("dd_dt", 277),
            ("dff_dt", 278),
            ("dfs_dt", 279),
            ("dfcaf_dt", 280),
            ("dfcas_dt", 281),
            ("djca_dt", 282),
            ("dffp_dt", 283),
            ("dfcafp_dt", 284),
            ("dnca_ss_dt", 285),
            ("dnca_i_dt", 286),
            ("dm_dt", 287),
            ("dh_dt", 288),
            ("dj_dt", 289),
            ("dhp_dt", 290),
            ("djp_dt", 291),
            ("dmL_dt", 292),
            ("dhL_dt", 293),
            ("dhLp_dt", 294),
            ("da_dt", 295),
            ("diF_dt", 296),
            ("diS_dt", 297),
            ("dap_dt", 298),
            ("diFp_dt", 299),
            ("diSp_dt", 300),
            ("dC3_dt", 301),
            ("dC2_dt", 302),
            ("dC1_dt", 303),
            ("dO_dt", 304),
            ("dI_dt", 305),
            ("dxs1_dt", 306),
            ("dxs2_dt", 307),
            ("dJrel_np_dt", 308),
            ("dJrel_p_dt", 309),
            ("dv_dt", 310),
            ("dnai_dt", 311),
            ("dnass_dt", 312),
            ("dki_dt", 313),
            ("dkss_dt", 314),
            ("dcli_dt", 315),
            ("dclss_dt", 316),
            ("dcai_dt", 317),
            ("dcass_dt", 318),
            ("dcansr_dt", 319),
            ("dcajsr_dt", 320),
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
    Compute the right hand side of the ToRORd_dyn_chloride_epi ODE
    """

    # Assign states
    assert len(states) == 45
    (
        CaMKt,
        d,
        ff,
        fs,
        fcaf,
        fcas,
        jca,
        ffp,
        fcafp,
        nca_ss,
        nca_i,
        m,
        h,
        j,
        hp,
        jp,
        mL,
        hL,
        hLp,
        a,
        iF,
        iS,
        ap,
        iFp,
        iSp,
        C3,
        C2,
        C1,
        O,
        I,
        xs1,
        xs2,
        Jrel_np,
        Jrel_p,
        v,
        nai,
        nass,
        ki,
        kss,
        cli,
        clss,
        cai,
        cass,
        cansr,
        cajsr,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    celltype = parameters[0]
    cao = parameters[1]
    clo = parameters[2]
    ko = parameters[3]
    nao = parameters[4]
    F = parameters[5]
    R = parameters[6]
    T = parameters[7]
    zca = parameters[8]
    zcl = parameters[9]
    zk = parameters[10]
    zna = parameters[11]
    L = parameters[12]
    rad = parameters[13]
    CaMKo = parameters[14]
    KmCaM = parameters[15]
    KmCaMK = parameters[16]
    aCaMK = parameters[17]
    bCaMK = parameters[18]
    PKNa = parameters[19]
    GNa = parameters[20]
    Gncx_b = parameters[21]
    INaCa_fractionSS = parameters[22]
    KmCaAct = parameters[23]
    kasymm = parameters[24]
    kcaoff = parameters[25]
    kcaon = parameters[26]
    kna1 = parameters[27]
    kna2 = parameters[28]
    kna3 = parameters[29]
    qca = parameters[30]
    qna = parameters[31]
    wca = parameters[32]
    wna = parameters[33]
    wnaca = parameters[34]
    H = parameters[35]
    Khp = parameters[36]
    Kki = parameters[37]
    Kko = parameters[38]
    Kmgatp = parameters[39]
    Knai0 = parameters[40]
    Knao0 = parameters[41]
    Knap = parameters[42]
    Kxkur = parameters[43]
    MgADP = parameters[44]
    MgATP = parameters[45]
    Pnak_b = parameters[46]
    delta = parameters[47]
    eP = parameters[48]
    k1m = parameters[49]
    k1p = parameters[50]
    k2m = parameters[51]
    k2p = parameters[52]
    k3m = parameters[53]
    k3p = parameters[54]
    k4m = parameters[55]
    k4p = parameters[56]
    PNab = parameters[57]
    GNaL_b = parameters[58]
    thL = parameters[59]
    GpCa = parameters[60]
    KmCap = parameters[61]
    tauCa = parameters[62]
    tauK = parameters[64]
    tauNa = parameters[65]
    Aff = parameters[66]
    ICaL_fractionSS = parameters[67]
    Kmn = parameters[68]
    PCa_b = parameters[69]
    dielConstant = parameters[70]
    k2n = parameters[71]
    offset = parameters[72]
    tjca = parameters[73]
    vShift = parameters[74]
    Jup_b = parameters[75]
    A_atp = parameters[76]
    K_atp = parameters[77]
    K_o_n = parameters[78]
    fkatp = parameters[79]
    gkatp = parameters[80]
    EKshift = parameters[81]
    Gto_b = parameters[82]
    GKr_b = parameters[83]
    alpha_1 = parameters[84]
    beta_1 = parameters[85]
    GKs_b = parameters[86]
    GK1_b = parameters[87]
    GKb_b = parameters[88]
    Fjunc = parameters[89]
    GClCa = parameters[90]
    GClb = parameters[91]
    KdClCa = parameters[92]
    PCab = parameters[93]
    Jrel_b = parameters[94]
    bt = parameters[95]
    cajsr_half = parameters[96]
    i_Stim_Amplitude = parameters[97]
    i_Stim_Period = parameters[99]
    i_Stim_PulseDuration = parameters[100]
    i_Stim_Start = parameters[101]
    BSLmax = parameters[102]
    BSRmax = parameters[103]
    KmBSL = parameters[104]
    KmBSR = parameters[105]
    cmdnmax_b = parameters[106]
    csqnmax = parameters[107]
    kmcmdn = parameters[108]
    kmcsqn = parameters[109]
    kmtrpn = parameters[110]
    trpnmax = parameters[111]

    # Init return args
    if values is None:
        values = np.zeros((45,), dtype=np.float_)
    else:
        assert isinstance(values, np.ndarray) and values.shape == (45,)

    # Expressions for the ToRORd dyn chloride component
    vfrt = F * v / (R * T)
    vffrt = (F * F) * v / (R * T)

    # Expressions for the Cell geometry component
    vcell = 3140.0 * L * (rad * rad)
    Ageo = 6.28 * (rad * rad) + 6.28 * L * rad
    Acap = 2 * Ageo
    vmyo = 0.68 * vcell
    vnsr = 0.0552 * vcell
    vjsr = 0.0048 * vcell
    vss = 0.02 * vcell

    # Expressions for the CaMK component
    CaMKb = CaMKo * (1 - CaMKt) / (1 + KmCaM / cass)
    CaMKa = CaMKb + CaMKt
    values[0] = -bCaMK * CaMKt + aCaMK * (CaMKb + CaMKt) * CaMKb

    # Expressions for the Reversal potentials component
    ENa = R * T * np.log(nao / nai) / (F * zna)
    EK = R * T * np.log(ko / ki) / (F * zk)
    EKs = R * T * np.log((ko + PKNa * nao) / (PKNa * nai + ki)) / (F * zk)
    ECl = R * T * np.log(clo / cli) / (F * zcl)
    EClss = R * T * np.log(clo / clss) / (F * zcl)

    # Expressions for the Ca component
    hca = np.exp(qca * vfrt)
    hna = np.exp(qna * vfrt)
    h1_i = 1 + (1 + hna) * nai / kna3
    h2_i = hna * nai / (kna3 * h1_i)
    h3_i = 1.0 / h1_i
    h4_i = 1 + (1 + nai / kna2) * nai / kna1
    h5_i = (nai * nai) / (kna1 * kna2 * h4_i)
    h6_i = 1.0 / h4_i
    h7_i = 1 + nao * (1 + 1.0 / hna) / kna3
    h8_i = nao / (kna3 * h7_i * hna)
    h9_i = 1.0 / h7_i
    h10_i = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    h11_i = (nao * nao) / (kna1 * kna2 * h10_i)
    h12_i = 1.0 / h10_i
    k1_i = cao * kcaon * h12_i
    k2_i = kcaoff
    k3p_i = wca * h9_i
    k3pp_i = wnaca * h8_i
    k3_i = k3p_i + k3pp_i
    k4p_i = wca * h3_i / hca
    k4pp_i = wnaca * h2_i
    k4_i = k4p_i + k4pp_i
    k5_i = kcaoff
    k6_i = kcaon * cai * h6_i
    k7_i = wna * h2_i * h5_i
    k8_i = wna * h11_i * h8_i
    x1_i = (k2_i + k3_i) * k5_i * k7_i + (k6_i + k7_i) * k2_i * k4_i
    x2_i = (k1_i + k8_i) * k4_i * k6_i + (k4_i + k5_i) * k1_i * k7_i
    x3_i = (k2_i + k3_i) * k6_i * k8_i + (k6_i + k7_i) * k1_i * k3_i
    x4_i = (k1_i + k8_i) * k3_i * k5_i + (k4_i + k5_i) * k2_i * k8_i
    E1_i = x1_i / (x1_i + x2_i + x3_i + x4_i)
    E2_i = x2_i / (x1_i + x2_i + x3_i + x4_i)
    E3_i = x3_i / (x1_i + x2_i + x3_i + x4_i)
    E4_i = x4_i / (x1_i + x2_i + x3_i + x4_i)
    allo_i = 1.0 / (1 + (KmCaAct * KmCaAct) / (cai * cai))
    JncxNa_i = E3_i * k4pp_i - E2_i * k3pp_i - 3 * E1_i * k8_i + 3 * E4_i * k7_i
    JncxCa_i = E2_i * k2_i - E1_i * k1_i
    Gncx = np.where(
        celltype == 1, 1.1 * Gncx_b, np.where(celltype == 2, 1.4 * Gncx_b, Gncx_b)
    )
    INaCa_i = (1 - INaCa_fractionSS) * (zca * JncxCa_i + zna * JncxNa_i) * Gncx * allo_i
    h1_ss = 1 + (1 + hna) * nass / kna3
    h2_ss = hna * nass / (kna3 * h1_ss)
    h3_ss = 1.0 / h1_ss
    h4_ss = 1 + (1 + nass / kna2) * nass / kna1
    h5_ss = (nass * nass) / (kna1 * kna2 * h4_ss)
    h6_ss = 1.0 / h4_ss
    h7_ss = 1 + nao * (1 + 1.0 / hna) / kna3
    h8_ss = nao / (kna3 * h7_ss * hna)
    h9_ss = 1.0 / h7_ss
    h10_ss = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    h11_ss = (nao * nao) / (kna1 * kna2 * h10_ss)
    h12_ss = 1.0 / h10_ss
    k1_ss = cao * kcaon * h12_ss
    k2_ss = kcaoff
    k3p_ss = wca * h9_ss
    k3pp_ss = wnaca * h8_ss
    k3_ss = k3p_ss + k3pp_ss
    k4p_ss = wca * h3_ss / hca
    k4pp_ss = wnaca * h2_ss
    k4_ss = k4p_ss + k4pp_ss
    k5_ss = kcaoff
    k6_ss = kcaon * cass * h6_ss
    k7_ss = wna * h2_ss * h5_ss
    k8_ss = wna * h11_ss * h8_ss
    x1_ss = (k2_ss + k3_ss) * k5_ss * k7_ss + (k6_ss + k7_ss) * k2_ss * k4_ss
    x2_ss = (k1_ss + k8_ss) * k4_ss * k6_ss + (k4_ss + k5_ss) * k1_ss * k7_ss
    x3_ss = (k2_ss + k3_ss) * k6_ss * k8_ss + (k6_ss + k7_ss) * k1_ss * k3_ss
    x4_ss = (k1_ss + k8_ss) * k3_ss * k5_ss + (k4_ss + k5_ss) * k2_ss * k8_ss
    E1_ss = x1_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E2_ss = x2_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E3_ss = x3_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E4_ss = x4_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    allo_ss = 1.0 / (1 + (KmCaAct * KmCaAct) / (cass * cass))
    JncxNa_ss = (
        E3_ss * k4pp_ss - E2_ss * k3pp_ss - 3 * E1_ss * k8_ss + 3 * E4_ss * k7_ss
    )
    JncxCa_ss = E2_ss * k2_ss - E1_ss * k1_ss
    INaCa_ss = INaCa_fractionSS * (zca * JncxCa_ss + zna * JncxNa_ss) * Gncx * allo_ss

    # Expressions for the K component
    Knai = Knai0 * np.exp(delta * vfrt / 3)
    Knao = Knao0 * np.exp((1 - delta) * vfrt / 3)
    P = eP / (1 + H / Khp + nai / Knap + ki / Kxkur)
    a1 = (
        k1p
        * (nai * nai * nai)
        / (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
            * (Knai * Knai * Knai)
        )
    )
    b1 = MgADP * k1m
    a2 = k2p
    b2 = (
        k2m
        * (nao * nao * nao)
        / (
            (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
            * (Knao * Knao * Knao)
        )
    )
    a3 = (
        k3p
        * (ko * ko)
        / (
            (Kko * Kko)
            * (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
        )
    )
    b3 = H * k3m * P / (1 + MgATP / Kmgatp)
    a4 = MgATP * k4p / (Kmgatp * (1 + MgATP / Kmgatp))
    b4 = (
        k4m
        * (ki * ki)
        / (
            (Kki * Kki)
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
        )
    )
    x1 = a1 * a2 * a4 + a1 * a2 * b3 + a2 * b3 * b4 + b2 * b3 * b4
    x2 = a1 * a2 * a3 + a2 * a3 * b4 + a3 * b1 * b4 + b1 * b2 * b4
    x3 = a2 * a3 * a4 + a3 * a4 * b1 + a4 * b1 * b2 + b1 * b2 * b3
    x4 = a1 * a3 * a4 + a1 * a4 * b2 + a1 * b2 * b3 + b2 * b3 * b4
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    JnakNa = -3 * E2 * b3 + 3 * E1 * a3
    JnakK = -2 * E3 * a1 + 2 * E4 * b1
    Pnak = np.where(
        celltype == 1, 0.9 * Pnak_b, np.where(celltype == 2, 0.7 * Pnak_b, Pnak_b)
    )
    INaK = (zk * JnakK + zna * JnakNa) * Pnak

    # Expressions for the b component
    INab = PNab * (-nao + np.exp(vfrt) * nai) * vffrt / (-1 + np.exp(vfrt))

    # Expressions for the IpCa component
    IpCa = GpCa * cai / (KmCap + cai)

    # Expressions for the Diff component
    JdiffNa = (-nai + nass) / tauNa
    JdiffK = (-ki + kss) / tauK
    Jdiff = (-cai + cass) / tauCa
    JdiffCl = (-cli + clss) / tauNa

    # Expressions for the Trans flux component
    Jtr = -cajsr / 60 + cansr / 60

    # Expressions for the ICaL component
    dss = np.where(v >= 31.4978, 1, 1.0763 * np.exp(-1.007 * np.exp(-0.0829 * v)))
    td = (
        0.6
        + offset
        + 1.0
        / (
            3.5254214873653824 * np.exp(0.09 * vShift + 0.09 * v)
            + 0.7408182206817179 * np.exp(-0.05 * vShift - 0.05 * v)
        )
    )
    values[1] = (-d + dss) / td
    fss = 1.0 / (1 + 199.86038496778565 * np.exp(0.27056277056277056 * v))
    tff = 7 + 1.0 / (0.0045 * np.exp(-2 - v / 10) + 0.0045 * np.exp(2 + v / 10))
    tfs = 1000 + 1.0 / (
        3.5e-05 * np.exp(-5 / 4 - v / 4) + 3.5e-05 * np.exp(5 / 6 + v / 6)
    )
    Afs = 1 - Aff
    values[2] = (-ff + fss) / tff
    values[3] = (-fs + fss) / tfs
    f = Aff * ff + Afs * fs
    fcass = fss
    tfcaf = 7 + 1.0 / (0.04 * np.exp(-4 / 7 + v / 7) + 0.04 * np.exp(4 / 7 - v / 7))
    tfcas = 100 + 1.0 / (0.00012 * np.exp(-v / 3) + 0.00012 * np.exp(v / 7))
    Afcaf = 0.3 + 0.6 / (1 + np.exp(-1 + v / 10))
    Afcas = 1 - Afcaf
    values[4] = (-fcaf + fcass) / tfcaf
    values[5] = (-fcas + fcass) / tfcas
    fca = Afcaf * fcaf + Afcas * fcas
    jcass = 1.0 / (1.0 + 649.7401897235336 * np.exp(0.35821750967187277 * v))
    values[6] = (-jca + jcass) / tjca
    tffp = 2.5 * tff
    values[7] = (-ffp + fss) / tffp
    fp = Aff * ffp + Afs * fs
    tfcafp = 2.5 * tfcaf
    values[8] = (-fcafp + fcass) / tfcafp
    fcap = Afcaf * fcafp + Afcas * fcas
    km2n = jca
    anca_ss = 1.0 / (np.power(1 + Kmn / cass, 4) + k2n / km2n)
    values[9] = k2n * anca_ss - km2n * nca_ss
    Io = 0.0005 * clo + 0.0005 * ko + 0.0005 * nao + 0.002 * cao
    Iss = 0.0005 * clss + 0.0005 * kss + 0.0005 * nass + 0.002 * cass
    constA = 1820000.0 * np.power(T * dielConstant, -1.5)
    gamma_cass = np.exp(-4 * (-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_cao = np.exp(-4 * (-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    gamma_nass = np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_nao = np.exp(-(-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    gamma_kss = np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_ko = np.exp(-(-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    PhiCaL_ss = (
        4
        * (-cao * gamma_cao + cass * np.exp(2 * vfrt) * gamma_cass)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    PhiCaNa_ss = (
        (-nao * gamma_nao + np.exp(vfrt) * gamma_nass * nass)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    PhiCaK_ss = (
        (-ko * gamma_ko + np.exp(vfrt) * gamma_kss * kss) * vffrt / (-1 + np.exp(vfrt))
    )
    PCa = np.where(
        celltype == 1, 1.2 * PCa_b, np.where(celltype == 2, 2 * PCa_b, PCa_b)
    )
    PCap = 1.1 * PCa
    PCaNa = 0.00125 * PCa
    PCaK = 0.0003574 * PCa
    PCaNap = 0.00125 * PCap
    PCaKp = 0.0003574 * PCap
    fICaLp = 1.0 / (1 + KmCaMK / CaMKa)
    ICaL_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCa * PhiCaL_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCap * PhiCaL_ss * d * fICaLp
    )
    ICaNa_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaNa * PhiCaNa_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaNap * PhiCaNa_ss * d * fICaLp
    )
    ICaK_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaK * PhiCaK_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaKp * PhiCaK_ss * d * fICaLp
    )
    anca_i = 1.0 / (np.power(1 + Kmn / cai, 4) + k2n / km2n)
    values[10] = k2n * anca_i - km2n * nca_i
    Ii = 0.0005 * cli + 0.0005 * ki + 0.0005 * nai + 0.002 * cai
    gamma_cai = np.exp(-4 * (-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    gamma_nai = np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    gamma_ki = np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    PhiCaL_i = (
        4
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    PhiCaNa_i = (
        (-nao * gamma_nao + np.exp(vfrt) * gamma_nai * nai)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    PhiCaK_i = (
        (-ko * gamma_ko + np.exp(vfrt) * gamma_ki * ki) * vffrt / (-1 + np.exp(vfrt))
    )
    ICaL_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCa * PhiCaL_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCap * PhiCaL_i * d * fICaLp
    )
    ICaNa_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaNa * PhiCaNa_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaNap * PhiCaNa_i * d * fICaLp
    )
    ICaK_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaK * PhiCaK_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaKp * PhiCaK_i * d * fICaLp
    )
    ICaL = ICaL_i + ICaL_ss
    ICaNa = ICaNa_i + ICaNa_ss
    ICaK = ICaK_i + ICaK_ss

    # Expressions for the SERCA component
    upScale = np.where(celltype == 1, 1.3, 1)
    Jupnp = 0.005425 * cai * upScale / (0.00092 + cai)
    Jupp = 0.01491875 * cai * upScale / (0.00075 + cai)
    fJupp = 1.0 / (1 + KmCaMK / CaMKa)
    Jleak = 0.0003255 * cansr
    Jup = Jup_b * (-Jleak + (1 - fJupp) * Jupnp + Jupp * fJupp)

    # Expressions for the I_katp component
    akik = np.power(ko / K_o_n, 0.24)
    bkik = 1.0 / (1 + (A_atp * A_atp) / (K_atp * K_atp))
    I_katp = fkatp * gkatp * (-EK + v) * akik * bkik

    # Expressions for the INa component
    mss = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
    )
    tm = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    values[11] = (-m + mss) / tm
    hss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
    )
    ah = np.where(
        v >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * v)
    )
    bh = np.where(
        v >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * v)),
        310000.0 * np.exp(0.3485 * v) + 2.7 * np.exp(0.079 * v),
    )
    th = 1.0 / (ah + bh)
    values[12] = (-h + hss) / th
    aj = np.where(
        v >= -40,
        0,
        (37.78 + v)
        * (-25428.0 * np.exp(0.2444 * v) - 6.948e-06 * np.exp(-0.04391 * v))
        / (1 + 50262745825.95399 * np.exp(0.311 * v)),
    )
    bj = np.where(
        v >= -40,
        0.6 * np.exp(0.057 * v) / (1 + 0.040762203978366204 * np.exp(-0.1 * v)),
        0.02424
        * np.exp(-0.01052 * v)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * v)),
    )
    jss = hss
    tj = 1.0 / (aj + bj)
    values[13] = (-j + jss) / tj
    hssp = 1.0 / (
        (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
        * (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
    )
    values[14] = (-hp + hssp) / th
    tjp = 1.46 * tj
    values[15] = (-jp + jss) / tjp
    fINap = 1.0 / (1 + KmCaMK / CaMKa)
    INa = GNa * (m * m * m) * (-ENa + v) * ((1 - fINap) * h * j + fINap * hp * jp)

    # Expressions for the L component
    mLss = 1.0 / (1 + 0.000291579585635531 * np.exp(-0.18996960486322187 * v))
    tmL = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    values[16] = (-mL + mLss) / tmL
    hLss = 1.0 / (1 + 120578.15595522427 * np.exp(0.13354700854700854 * v))
    values[17] = (-hL + hLss) / thL
    hLssp = 1.0 / (1 + 275969.2903869871 * np.exp(0.13354700854700854 * v))
    thLp = 3 * thL
    values[18] = (-hLp + hLssp) / thLp
    GNaL = np.where(celltype == 1, 0.6 * GNaL_b, GNaL_b)
    fINaLp = 1.0 / (1 + KmCaMK / CaMKa)
    INaL = (-ENa + v) * ((1 - fINaLp) * hL + fINaLp * hLp) * GNaL * mL

    # Expressions for the Ito component
    ass = 1.0 / (
        1
        + 2.6316508161673635
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    ta = 1.0515 / (
        1.0
        / (
            1.2089
            + 2.2621017070578837
            * np.exp(-0.03403513787634354 * EKshift - 0.03403513787634354 * v)
        )
        + 3.5
        / (
            1
            + 30.069572727397507
            * np.exp(0.03403513787634354 * EKshift + 0.03403513787634354 * v)
        )
    )
    values[19] = (-a + ass) / ta
    iss = 1.0 / (
        1
        + 2194.970764538301
        * np.exp(0.17510068289266328 * EKshift + 0.17510068289266328 * v)
    )
    delta_epi = np.where(
        celltype == 1, 1 - 0.95 / (1 + np.exp(14 + EKshift / 5 + v / 5)), 1
    )
    tiF_b = 4.562 + 1.0 / (
        0.3933 * np.exp(-1 - EKshift / 100 - v / 100)
        + 1.6300896349780942
        * np.exp(0.06027727546714889 * EKshift + 0.06027727546714889 * v)
    )
    tiS_b = 23.62 + 1.0 / (
        0.00027617763953377436
        * np.exp(-0.01693480101608806 * EKshift - 0.01693480101608806 * v)
        + 0.024208962804604526
        * np.exp(0.12377769525931426 * EKshift + 0.12377769525931426 * v)
    )
    tiF = delta_epi * tiF_b
    tiS = delta_epi * tiS_b
    AiF = 1.0 / (
        1
        + 0.24348537187522867
        * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
    )
    AiS = 1 - AiF
    values[20] = (-iF + iss) / tiF
    values[21] = (-iS + iss) / tiS
    i = AiF * iF + AiS * iS
    assp = 1.0 / (
        1
        + 5.167428462230666
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    values[22] = (-ap + assp) / ta
    dti_develop = 1.354 + 0.0001 / (
        2.6591269045230603e-05
        * np.exp(0.06293266205160478 * EKshift + 0.06293266205160478 * v)
        + 4.5541779737128264e24
        * np.exp(-4.642525533890436 * EKshift - 4.642525533890436 * v)
    )
    dti_recover = 1 - 0.5 / (1 + 33.11545195869231 * np.exp(0.05 * EKshift + 0.05 * v))
    tiFp = dti_develop * dti_recover * tiF
    tiSp = dti_develop * dti_recover * tiS
    values[23] = (-iFp + iss) / tiFp
    values[24] = (-iSp + iss) / tiSp
    ip = AiF * iFp + AiS * iSp
    Gto = np.where(celltype == 1, 2 * Gto_b, np.where(celltype == 2, 2 * Gto_b, Gto_b))
    fItop = 1.0 / (1 + KmCaMK / CaMKa)
    Ito = (-EK + v) * ((1 - fItop) * a * i + ap * fItop * ip) * Gto

    # Expressions for the IKr component
    alpha = 0.1161 * np.exp(0.299 * vfrt)
    beta = 0.2442 * np.exp(-1.604 * vfrt)
    alpha_2 = 0.0578 * np.exp(0.971 * vfrt)
    beta_2 = 0.000349 * np.exp(-1.062 * vfrt)
    alpha_i = 0.2533 * np.exp(0.5953 * vfrt)
    beta_i = 0.06525 * np.exp(-0.8209 * vfrt)
    alpha_C2ToI = 5.2e-05 * np.exp(1.525 * vfrt)
    beta_ItoC2 = alpha_C2ToI * beta_2 * beta_i / (alpha_2 * alpha_i)
    values[25] = C2 * beta - C3 * alpha
    values[26] = beta_1 * C1 + C3 * alpha - (alpha_1 + beta) * C2
    values[27] = (
        alpha_1 * C2
        + I * beta_ItoC2
        + O * beta_2
        - (beta_1 + alpha_2 + alpha_C2ToI) * C1
    )
    values[28] = C1 * alpha_2 + I * beta_i - (alpha_i + beta_2) * O
    values[29] = C1 * alpha_C2ToI + O * alpha_i - (beta_ItoC2 + beta_i) * I
    GKr = np.where(
        celltype == 1, 1.3 * GKr_b, np.where(celltype == 2, 0.8 * GKr_b, GKr_b)
    )
    IKr = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GKr * O / 5

    # Expressions for the IKs component
    xs1ss = 1.0 / (1 + 0.27288596035656526 * np.exp(-0.11195700850873264 * v))
    txs1 = 817.3 + 1.0 / (
        0.003504067763074858 * np.exp(0.056179775280898875 * v)
        + 0.001292 * np.exp(-21 / 23 - v / 230)
    )
    values[30] = (-xs1 + xs1ss) / txs1
    xs2ss = xs1ss
    txs2 = 1.0 / (
        0.0022561357010639103 * np.exp(-v / 31) + 0.01 * np.exp(-5 / 2 + v / 20)
    )
    values[31] = (-xs2 + xs2ss) / txs2
    KsCa = 1 + 0.6 / (1 + 6.481821026062645e-07 * np.power(1.0 / cai, 1.4))
    GKs = np.where(celltype == 1, 1.4 * GKs_b, GKs_b)
    IKs = (-EKs + v) * GKs * KsCa * xs1 * xs2

    # Expressions for the IK1 component
    aK1 = 4.094 / (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
    bK1 = (
        12.621629724407278 * np.exp(0.0674 * v - 0.0674 * EK)
        + 1.1196358381249121e-16 * np.exp(0.0618 * v - 0.0618 * EK)
    ) / (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
    K1ss = aK1 / (aK1 + bK1)
    GK1 = np.where(
        celltype == 1, 1.2 * GK1_b, np.where(celltype == 2, 1.3 * GK1_b, GK1_b)
    )
    IK1 = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GK1 * K1ss / 5

    # Expressions for the IKb component
    xkb = 1.0 / (1 + 1.57503502085457 * np.exp(-0.04168907454423419 * v))
    GKb = np.where(celltype == 1, 0.6 * GKb_b, GKb_b)
    IKb = (-EK + v) * GKb * xkb

    # Expressions for the ICl component
    IClCa_junc = Fjunc * GClCa * (-EClss + v) / (1 + KdClCa / cass)
    IClCa_sl = GClCa * (1 - Fjunc) * (-ECl + v) / (1 + KdClCa / cai)
    IClCa = IClCa_junc + IClCa_sl
    IClb = GClb * (-ECl + v)

    # Expressions for the ICab component
    ICab = (
        4
        * PCab
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )

    # Expressions for the Ryr component
    a_rel = 0.5 * bt
    Jrel_inf_b = -ICaL_ss * a_rel / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    Jrel_inf = np.where(celltype == 2, 1.7 * Jrel_inf_b, Jrel_inf_b)
    tau_rel_b = bt / (1 + 0.0123 / cajsr)
    tau_rel = np.where(tau_rel_b < 0.001, 0.001, tau_rel_b)
    values[32] = (-Jrel_np + Jrel_inf) / tau_rel
    btp = 1.25 * bt
    a_relp = 0.5 * btp
    Jrel_infp_b = -ICaL_ss * a_relp / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    Jrel_infp = np.where(celltype == 2, 1.7 * Jrel_infp_b, Jrel_infp_b)
    tau_relp_b = btp / (1 + 0.0123 / cajsr)
    tau_relp = np.where(tau_relp_b < 0.001, 0.001, tau_relp_b)
    values[33] = (-Jrel_p + Jrel_infp) / tau_relp
    fJrelp = 1.0 / (1 + KmCaMK / CaMKa)
    Jrel = Jrel_b * ((1 - fJrelp) * Jrel_np + Jrel_p * fJrelp)

    # Expressions for the Membrane component
    Istim = np.where(
        np.logical_and(
            t >= i_Stim_Start,
            t
            - i_Stim_Start
            - i_Stim_Period * np.floor((t - i_Stim_Start) / i_Stim_Period)
            <= i_Stim_PulseDuration,
        ),
        i_Stim_Amplitude,
        0,
    )
    values[34] = (
        -ICaK
        - ICaL
        - ICaNa
        - ICab
        - IClCa
        - IClb
        - IK1
        - IKb
        - IKr
        - IKs
        - INa
        - INaCa_i
        - INaCa_ss
        - INaK
        - INaL
        - INab
        - I_katp
        - IpCa
        - Istim
        - Ito
    )

    # Expressions for the Intracellular ions component
    cmdnmax = np.where(celltype == 1, 1.3 * cmdnmax_b, cmdnmax_b)
    values[35] = JdiffNa * vss / vmyo + (
        -ICaNa_i - INa - INaL - INab - 3 * INaCa_i - 3 * INaK
    ) * Acap / (F * vmyo)
    values[36] = -JdiffNa + (-ICaNa_ss - 3 * INaCa_ss) * Acap / (F * vss)
    values[37] = JdiffK * vss / vmyo + (
        -ICaK_i - IK1 - IKb - IKr - IKs - I_katp - Istim - Ito + 2 * INaK
    ) * Acap / (F * vmyo)
    values[38] = -JdiffK - Acap * ICaK_ss / (F * vss)
    values[39] = JdiffCl * vss / vmyo + (IClCa_sl + IClb) * Acap / (F * vmyo)
    values[40] = -JdiffCl + Acap * IClCa_junc / (F * vss)
    Bcai = 1.0 / (
        1
        + kmcmdn * cmdnmax / ((kmcmdn + cai) * (kmcmdn + cai))
        + kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai))
    )
    values[41] = (
        Jdiff * vss / vmyo
        - Jup * vnsr / vmyo
        + (-ICaL_i - ICab - IpCa + 2 * INaCa_i) * Acap / (2 * F * vmyo)
    ) * Bcai
    Bcass = 1.0 / (
        1
        + BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass))
        + BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass))
    )
    values[42] = (
        -Jdiff + Jrel * vjsr / vss + (-ICaL_ss + 2 * INaCa_ss) * Acap / (2 * F * vss)
    ) * Bcass
    values[43] = -Jtr * vjsr / vnsr + Jup
    Bcajsr = 1.0 / (1 + csqnmax * kmcsqn / ((kmcsqn + cajsr) * (kmcsqn + cajsr)))
    values[44] = (-Jrel + Jtr) * Bcajsr

    # Return results
    return values


def monitor(states, t, parameters, monitored=None):
    """
    Computes monitored expressions of the ToRORd_dyn_chloride_epi ODE
    """

    # Assign states
    assert len(states) == 45
    (
        CaMKt,
        d,
        ff,
        fs,
        fcaf,
        fcas,
        jca,
        ffp,
        fcafp,
        nca_ss,
        nca_i,
        m,
        h,
        j,
        hp,
        jp,
        mL,
        hL,
        hLp,
        a,
        iF,
        iS,
        ap,
        iFp,
        iSp,
        C3,
        C2,
        C1,
        O,
        I,
        xs1,
        xs2,
        Jrel_np,
        Jrel_p,
        v,
        nai,
        nass,
        ki,
        kss,
        cli,
        clss,
        cai,
        cass,
        cansr,
        cajsr,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    celltype = parameters[0]
    cao = parameters[1]
    clo = parameters[2]
    ko = parameters[3]
    nao = parameters[4]
    F = parameters[5]
    R = parameters[6]
    T = parameters[7]
    zca = parameters[8]
    zcl = parameters[9]
    zk = parameters[10]
    zna = parameters[11]
    L = parameters[12]
    rad = parameters[13]
    CaMKo = parameters[14]
    KmCaM = parameters[15]
    KmCaMK = parameters[16]
    aCaMK = parameters[17]
    bCaMK = parameters[18]
    PKNa = parameters[19]
    GNa = parameters[20]
    Gncx_b = parameters[21]
    INaCa_fractionSS = parameters[22]
    KmCaAct = parameters[23]
    kasymm = parameters[24]
    kcaoff = parameters[25]
    kcaon = parameters[26]
    kna1 = parameters[27]
    kna2 = parameters[28]
    kna3 = parameters[29]
    qca = parameters[30]
    qna = parameters[31]
    wca = parameters[32]
    wna = parameters[33]
    wnaca = parameters[34]
    H = parameters[35]
    Khp = parameters[36]
    Kki = parameters[37]
    Kko = parameters[38]
    Kmgatp = parameters[39]
    Knai0 = parameters[40]
    Knao0 = parameters[41]
    Knap = parameters[42]
    Kxkur = parameters[43]
    MgADP = parameters[44]
    MgATP = parameters[45]
    Pnak_b = parameters[46]
    delta = parameters[47]
    eP = parameters[48]
    k1m = parameters[49]
    k1p = parameters[50]
    k2m = parameters[51]
    k2p = parameters[52]
    k3m = parameters[53]
    k3p = parameters[54]
    k4m = parameters[55]
    k4p = parameters[56]
    PNab = parameters[57]
    GNaL_b = parameters[58]
    thL = parameters[59]
    GpCa = parameters[60]
    KmCap = parameters[61]
    tauCa = parameters[62]
    tauK = parameters[64]
    tauNa = parameters[65]
    Aff = parameters[66]
    ICaL_fractionSS = parameters[67]
    Kmn = parameters[68]
    PCa_b = parameters[69]
    dielConstant = parameters[70]
    k2n = parameters[71]
    offset = parameters[72]
    tjca = parameters[73]
    vShift = parameters[74]
    Jup_b = parameters[75]
    A_atp = parameters[76]
    K_atp = parameters[77]
    K_o_n = parameters[78]
    fkatp = parameters[79]
    gkatp = parameters[80]
    EKshift = parameters[81]
    Gto_b = parameters[82]
    GKr_b = parameters[83]
    alpha_1 = parameters[84]
    beta_1 = parameters[85]
    GKs_b = parameters[86]
    GK1_b = parameters[87]
    GKb_b = parameters[88]
    Fjunc = parameters[89]
    GClCa = parameters[90]
    GClb = parameters[91]
    KdClCa = parameters[92]
    PCab = parameters[93]
    Jrel_b = parameters[94]
    bt = parameters[95]
    cajsr_half = parameters[96]
    i_Stim_Amplitude = parameters[97]
    i_Stim_Period = parameters[99]
    i_Stim_PulseDuration = parameters[100]
    i_Stim_Start = parameters[101]
    BSLmax = parameters[102]
    BSRmax = parameters[103]
    KmBSL = parameters[104]
    KmBSR = parameters[105]
    cmdnmax_b = parameters[106]
    csqnmax = parameters[107]
    kmcmdn = parameters[108]
    kmcsqn = parameters[109]
    kmtrpn = parameters[110]
    trpnmax = parameters[111]

    # Init return args
    if monitored is None:
        monitored = np.zeros((321,), dtype=np.float_)
    else:
        assert isinstance(monitored, np.ndarray) and monitored.shape == (321,)

    # Expressions for the ToRORd dyn chloride component
    monitored[273] = F * v / (R * T)
    monitored[274] = (F * F) * v / (R * T)

    # Expressions for the Cell geometry component
    monitored[0] = 3140.0 * L * (rad * rad)
    monitored[1] = 6.28 * (rad * rad) + 6.28 * L * rad
    monitored[2] = 2 * monitored[1]
    monitored[3] = 0.68 * monitored[0]
    monitored[4] = 0.0552 * monitored[0]
    monitored[5] = 0.0048 * monitored[0]
    monitored[6] = 0.02 * monitored[0]

    # Expressions for the CaMK component
    monitored[7] = CaMKo * (1 - CaMKt) / (1 + KmCaM / cass)
    monitored[8] = CaMKt + monitored[7]
    monitored[276] = -bCaMK * CaMKt + aCaMK * (CaMKt + monitored[7]) * monitored[7]

    # Expressions for the Reversal potentials component
    monitored[9] = R * T * np.log(nao / nai) / (F * zna)
    monitored[10] = R * T * np.log(ko / ki) / (F * zk)
    monitored[11] = R * T * np.log((ko + PKNa * nao) / (PKNa * nai + ki)) / (F * zk)
    monitored[12] = R * T * np.log(clo / cli) / (F * zcl)
    monitored[13] = R * T * np.log(clo / clss) / (F * zcl)

    # Expressions for the Ca component
    monitored[28] = np.exp(qca * monitored[273])
    monitored[29] = np.exp(qna * monitored[273])
    monitored[30] = 1 + (1 + monitored[29]) * nai / kna3
    monitored[31] = monitored[29] * nai / (kna3 * monitored[30])
    monitored[32] = 1.0 / monitored[30]
    monitored[33] = 1 + (1 + nai / kna2) * nai / kna1
    monitored[34] = (nai * nai) / (kna1 * kna2 * monitored[33])
    monitored[35] = 1.0 / monitored[33]
    monitored[36] = 1 + nao * (1 + 1.0 / monitored[29]) / kna3
    monitored[37] = nao / (kna3 * monitored[29] * monitored[36])
    monitored[38] = 1.0 / monitored[36]
    monitored[39] = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    monitored[40] = (nao * nao) / (kna1 * kna2 * monitored[39])
    monitored[41] = 1.0 / monitored[39]
    monitored[42] = cao * kcaon * monitored[41]
    monitored[43] = kcaoff
    monitored[44] = wca * monitored[38]
    monitored[45] = wnaca * monitored[37]
    monitored[46] = monitored[44] + monitored[45]
    monitored[47] = wca * monitored[32] / monitored[28]
    monitored[48] = wnaca * monitored[31]
    monitored[49] = monitored[47] + monitored[48]
    monitored[50] = kcaoff
    monitored[51] = kcaon * cai * monitored[35]
    monitored[52] = wna * monitored[31] * monitored[34]
    monitored[53] = wna * monitored[37] * monitored[40]
    monitored[54] = (monitored[43] + monitored[46]) * monitored[50] * monitored[52] + (
        monitored[51] + monitored[52]
    ) * monitored[43] * monitored[49]
    monitored[55] = (monitored[42] + monitored[53]) * monitored[49] * monitored[51] + (
        monitored[49] + monitored[50]
    ) * monitored[42] * monitored[52]
    monitored[56] = (monitored[43] + monitored[46]) * monitored[51] * monitored[53] + (
        monitored[51] + monitored[52]
    ) * monitored[42] * monitored[46]
    monitored[57] = (monitored[42] + monitored[53]) * monitored[46] * monitored[50] + (
        monitored[49] + monitored[50]
    ) * monitored[43] * monitored[53]
    monitored[58] = monitored[54] / (
        monitored[54] + monitored[55] + monitored[56] + monitored[57]
    )
    monitored[59] = monitored[55] / (
        monitored[54] + monitored[55] + monitored[56] + monitored[57]
    )
    monitored[60] = monitored[56] / (
        monitored[54] + monitored[55] + monitored[56] + monitored[57]
    )
    monitored[61] = monitored[57] / (
        monitored[54] + monitored[55] + monitored[56] + monitored[57]
    )
    monitored[62] = 1.0 / (1 + (KmCaAct * KmCaAct) / (cai * cai))
    monitored[63] = (
        monitored[48] * monitored[60]
        - monitored[45] * monitored[59]
        - 3 * monitored[53] * monitored[58]
        + 3 * monitored[52] * monitored[61]
    )
    monitored[64] = monitored[43] * monitored[59] - monitored[42] * monitored[58]
    monitored[65] = np.where(
        celltype == 1, 1.1 * Gncx_b, np.where(celltype == 2, 1.4 * Gncx_b, Gncx_b)
    )
    monitored[66] = (
        (1 - INaCa_fractionSS)
        * (zca * monitored[64] + zna * monitored[63])
        * monitored[62]
        * monitored[65]
    )
    monitored[67] = 1 + (1 + monitored[29]) * nass / kna3
    monitored[68] = monitored[29] * nass / (kna3 * monitored[67])
    monitored[69] = 1.0 / monitored[67]
    monitored[70] = 1 + (1 + nass / kna2) * nass / kna1
    monitored[71] = (nass * nass) / (kna1 * kna2 * monitored[70])
    monitored[72] = 1.0 / monitored[70]
    monitored[73] = 1 + nao * (1 + 1.0 / monitored[29]) / kna3
    monitored[74] = nao / (kna3 * monitored[29] * monitored[73])
    monitored[75] = 1.0 / monitored[73]
    monitored[76] = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    monitored[77] = (nao * nao) / (kna1 * kna2 * monitored[76])
    monitored[78] = 1.0 / monitored[76]
    monitored[79] = cao * kcaon * monitored[78]
    monitored[80] = kcaoff
    monitored[81] = wca * monitored[75]
    monitored[82] = wnaca * monitored[74]
    monitored[83] = monitored[81] + monitored[82]
    monitored[84] = wca * monitored[69] / monitored[28]
    monitored[85] = wnaca * monitored[68]
    monitored[86] = monitored[84] + monitored[85]
    monitored[87] = kcaoff
    monitored[88] = kcaon * cass * monitored[72]
    monitored[89] = wna * monitored[68] * monitored[71]
    monitored[90] = wna * monitored[74] * monitored[77]
    monitored[91] = (monitored[80] + monitored[83]) * monitored[87] * monitored[89] + (
        monitored[88] + monitored[89]
    ) * monitored[80] * monitored[86]
    monitored[92] = (monitored[79] + monitored[90]) * monitored[86] * monitored[88] + (
        monitored[86] + monitored[87]
    ) * monitored[79] * monitored[89]
    monitored[93] = (monitored[80] + monitored[83]) * monitored[88] * monitored[90] + (
        monitored[88] + monitored[89]
    ) * monitored[79] * monitored[83]
    monitored[94] = (monitored[79] + monitored[90]) * monitored[83] * monitored[87] + (
        monitored[86] + monitored[87]
    ) * monitored[80] * monitored[90]
    monitored[95] = monitored[91] / (
        monitored[91] + monitored[92] + monitored[93] + monitored[94]
    )
    monitored[96] = monitored[92] / (
        monitored[91] + monitored[92] + monitored[93] + monitored[94]
    )
    monitored[97] = monitored[93] / (
        monitored[91] + monitored[92] + monitored[93] + monitored[94]
    )
    monitored[98] = monitored[94] / (
        monitored[91] + monitored[92] + monitored[93] + monitored[94]
    )
    monitored[99] = 1.0 / (1 + (KmCaAct * KmCaAct) / (cass * cass))
    monitored[100] = (
        monitored[85] * monitored[97]
        - monitored[82] * monitored[96]
        - 3 * monitored[90] * monitored[95]
        + 3 * monitored[89] * monitored[98]
    )
    monitored[101] = monitored[80] * monitored[96] - monitored[79] * monitored[95]
    monitored[102] = (
        INaCa_fractionSS
        * (zca * monitored[101] + zna * monitored[100])
        * monitored[65]
        * monitored[99]
    )

    # Expressions for the K component
    monitored[103] = Knai0 * np.exp(delta * monitored[273] / 3)
    monitored[104] = Knao0 * np.exp((1 - delta) * monitored[273] / 3)
    monitored[105] = eP / (1 + H / Khp + nai / Knap + ki / Kxkur)
    monitored[106] = (
        k1p
        * (nai * nai * nai)
        / (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + (
                    (1 + nai / monitored[103])
                    * (1 + nai / monitored[103])
                    * (1 + nai / monitored[103])
                )
            )
            * (monitored[103] * monitored[103] * monitored[103])
        )
    )
    monitored[107] = MgADP * k1m
    monitored[108] = k2p
    monitored[109] = (
        k2m
        * (nao * nao * nao)
        / (
            (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + (
                    (1 + nao / monitored[104])
                    * (1 + nao / monitored[104])
                    * (1 + nao / monitored[104])
                )
            )
            * (monitored[104] * monitored[104] * monitored[104])
        )
    )
    monitored[110] = (
        k3p
        * (ko * ko)
        / (
            (Kko * Kko)
            * (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + (
                    (1 + nao / monitored[104])
                    * (1 + nao / monitored[104])
                    * (1 + nao / monitored[104])
                )
            )
        )
    )
    monitored[111] = H * k3m * monitored[105] / (1 + MgATP / Kmgatp)
    monitored[112] = MgATP * k4p / (Kmgatp * (1 + MgATP / Kmgatp))
    monitored[113] = (
        k4m
        * (ki * ki)
        / (
            (Kki * Kki)
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + (
                    (1 + nai / monitored[103])
                    * (1 + nai / monitored[103])
                    * (1 + nai / monitored[103])
                )
            )
        )
    )
    monitored[114] = (
        monitored[106] * monitored[108] * monitored[111]
        + monitored[106] * monitored[108] * monitored[112]
        + monitored[108] * monitored[111] * monitored[113]
        + monitored[109] * monitored[111] * monitored[113]
    )
    monitored[115] = (
        monitored[106] * monitored[108] * monitored[110]
        + monitored[107] * monitored[109] * monitored[113]
        + monitored[107] * monitored[110] * monitored[113]
        + monitored[108] * monitored[110] * monitored[113]
    )
    monitored[116] = (
        monitored[107] * monitored[109] * monitored[111]
        + monitored[107] * monitored[109] * monitored[112]
        + monitored[107] * monitored[110] * monitored[112]
        + monitored[108] * monitored[110] * monitored[112]
    )
    monitored[117] = (
        monitored[106] * monitored[109] * monitored[111]
        + monitored[106] * monitored[109] * monitored[112]
        + monitored[106] * monitored[110] * monitored[112]
        + monitored[109] * monitored[111] * monitored[113]
    )
    monitored[118] = monitored[114] / (
        monitored[114] + monitored[115] + monitored[116] + monitored[117]
    )
    monitored[119] = monitored[115] / (
        monitored[114] + monitored[115] + monitored[116] + monitored[117]
    )
    monitored[120] = monitored[116] / (
        monitored[114] + monitored[115] + monitored[116] + monitored[117]
    )
    monitored[121] = monitored[117] / (
        monitored[114] + monitored[115] + monitored[116] + monitored[117]
    )
    monitored[122] = (
        -3 * monitored[111] * monitored[119] + 3 * monitored[110] * monitored[118]
    )
    monitored[123] = (
        -2 * monitored[106] * monitored[120] + 2 * monitored[107] * monitored[121]
    )
    monitored[124] = np.where(
        celltype == 1, 0.9 * Pnak_b, np.where(celltype == 2, 0.7 * Pnak_b, Pnak_b)
    )
    monitored[125] = (zk * monitored[123] + zna * monitored[122]) * monitored[124]

    # Expressions for the b component
    monitored[126] = (
        PNab
        * (-nao + np.exp(monitored[273]) * nai)
        * monitored[274]
        / (-1 + np.exp(monitored[273]))
    )

    # Expressions for the IpCa component
    monitored[135] = GpCa * cai / (KmCap + cai)

    # Expressions for the Diff component
    monitored[136] = (-nai + nass) / tauNa
    monitored[137] = (-ki + kss) / tauK
    monitored[138] = (-cai + cass) / tauCa
    monitored[139] = (-cli + clss) / tauNa

    # Expressions for the Trans flux component
    monitored[275] = -cajsr / 60 + cansr / 60

    # Expressions for the ICaL component
    monitored[140] = np.where(
        v >= 31.4978, 1, 1.0763 * np.exp(-1.007 * np.exp(-0.0829 * v))
    )
    monitored[141] = (
        0.6
        + offset
        + 1.0
        / (
            3.5254214873653824 * np.exp(0.09 * vShift + 0.09 * v)
            + 0.7408182206817179 * np.exp(-0.05 * vShift - 0.05 * v)
        )
    )
    monitored[277] = (-d + monitored[140]) / monitored[141]
    monitored[142] = 1.0 / (1 + 199.86038496778565 * np.exp(0.27056277056277056 * v))
    monitored[143] = 7 + 1.0 / (
        0.0045 * np.exp(-2 - v / 10) + 0.0045 * np.exp(2 + v / 10)
    )
    monitored[144] = 1000 + 1.0 / (
        3.5e-05 * np.exp(-5 / 4 - v / 4) + 3.5e-05 * np.exp(5 / 6 + v / 6)
    )
    monitored[145] = 1 - Aff
    monitored[278] = (-ff + monitored[142]) / monitored[143]
    monitored[279] = (-fs + monitored[142]) / monitored[144]
    monitored[146] = Aff * ff + fs * monitored[145]
    monitored[147] = monitored[142]
    monitored[148] = 7 + 1.0 / (
        0.04 * np.exp(-4 / 7 + v / 7) + 0.04 * np.exp(4 / 7 - v / 7)
    )
    monitored[149] = 100 + 1.0 / (0.00012 * np.exp(-v / 3) + 0.00012 * np.exp(v / 7))
    monitored[150] = 0.3 + 0.6 / (1 + np.exp(-1 + v / 10))
    monitored[151] = 1 - monitored[150]
    monitored[280] = (-fcaf + monitored[147]) / monitored[148]
    monitored[281] = (-fcas + monitored[147]) / monitored[149]
    monitored[152] = fcaf * monitored[150] + fcas * monitored[151]
    monitored[153] = 1.0 / (1.0 + 649.7401897235336 * np.exp(0.35821750967187277 * v))
    monitored[282] = (-jca + monitored[153]) / tjca
    monitored[154] = 2.5 * monitored[143]
    monitored[283] = (-ffp + monitored[142]) / monitored[154]
    monitored[155] = Aff * ffp + fs * monitored[145]
    monitored[156] = 2.5 * monitored[148]
    monitored[284] = (-fcafp + monitored[147]) / monitored[156]
    monitored[157] = fcafp * monitored[150] + fcas * monitored[151]
    monitored[158] = jca
    monitored[159] = 1.0 / (np.power(1 + Kmn / cass, 4) + k2n / monitored[158])
    monitored[285] = k2n * monitored[159] - monitored[158] * nca_ss
    monitored[160] = 0.0005 * clo + 0.0005 * ko + 0.0005 * nao + 0.002 * cao
    monitored[161] = 0.0005 * clss + 0.0005 * kss + 0.0005 * nass + 0.002 * cass
    monitored[162] = 1820000.0 * np.power(T * dielConstant, -1.5)
    monitored[163] = np.exp(
        -4
        * (
            -0.3 * monitored[161]
            + np.sqrt(monitored[161]) / (1 + np.sqrt(monitored[161]))
        )
        * monitored[162]
    )
    monitored[164] = np.exp(
        -4
        * (
            -0.3 * monitored[160]
            + np.sqrt(monitored[160]) / (1 + np.sqrt(monitored[160]))
        )
        * monitored[162]
    )
    monitored[165] = np.exp(
        -(
            -0.3 * monitored[161]
            + np.sqrt(monitored[161]) / (1 + np.sqrt(monitored[161]))
        )
        * monitored[162]
    )
    monitored[166] = np.exp(
        -(
            -0.3 * monitored[160]
            + np.sqrt(monitored[160]) / (1 + np.sqrt(monitored[160]))
        )
        * monitored[162]
    )
    monitored[167] = np.exp(
        -(
            -0.3 * monitored[161]
            + np.sqrt(monitored[161]) / (1 + np.sqrt(monitored[161]))
        )
        * monitored[162]
    )
    monitored[168] = np.exp(
        -(
            -0.3 * monitored[160]
            + np.sqrt(monitored[160]) / (1 + np.sqrt(monitored[160]))
        )
        * monitored[162]
    )
    monitored[169] = (
        4
        * (-cao * monitored[164] + cass * np.exp(2 * monitored[273]) * monitored[163])
        * monitored[274]
        / (-1 + np.exp(2 * monitored[273]))
    )
    monitored[170] = (
        (-nao * monitored[166] + np.exp(monitored[273]) * monitored[165] * nass)
        * monitored[274]
        / (-1 + np.exp(monitored[273]))
    )
    monitored[171] = (
        (-ko * monitored[168] + np.exp(monitored[273]) * kss * monitored[167])
        * monitored[274]
        / (-1 + np.exp(monitored[273]))
    )
    monitored[172] = np.where(
        celltype == 1, 1.2 * PCa_b, np.where(celltype == 2, 2 * PCa_b, PCa_b)
    )
    monitored[173] = 1.1 * monitored[172]
    monitored[174] = 0.00125 * monitored[172]
    monitored[175] = 0.0003574 * monitored[172]
    monitored[176] = 0.00125 * monitored[173]
    monitored[177] = 0.0003574 * monitored[173]
    monitored[178] = 1.0 / (1 + KmCaMK / monitored[8])
    monitored[179] = ICaL_fractionSS * (
        (1 - monitored[178])
        * ((1 - nca_ss) * monitored[146] + jca * monitored[152] * nca_ss)
        * d
        * monitored[169]
        * monitored[172]
        + ((1 - nca_ss) * monitored[155] + jca * monitored[157] * nca_ss)
        * d
        * monitored[169]
        * monitored[173]
        * monitored[178]
    )
    monitored[180] = ICaL_fractionSS * (
        (1 - monitored[178])
        * ((1 - nca_ss) * monitored[146] + jca * monitored[152] * nca_ss)
        * d
        * monitored[170]
        * monitored[174]
        + ((1 - nca_ss) * monitored[155] + jca * monitored[157] * nca_ss)
        * d
        * monitored[170]
        * monitored[176]
        * monitored[178]
    )
    monitored[181] = ICaL_fractionSS * (
        (1 - monitored[178])
        * ((1 - nca_ss) * monitored[146] + jca * monitored[152] * nca_ss)
        * d
        * monitored[171]
        * monitored[175]
        + ((1 - nca_ss) * monitored[155] + jca * monitored[157] * nca_ss)
        * d
        * monitored[171]
        * monitored[177]
        * monitored[178]
    )
    monitored[182] = 1.0 / (np.power(1 + Kmn / cai, 4) + k2n / monitored[158])
    monitored[286] = k2n * monitored[182] - monitored[158] * nca_i
    monitored[183] = 0.0005 * cli + 0.0005 * ki + 0.0005 * nai + 0.002 * cai
    monitored[184] = np.exp(
        -4
        * (
            -0.3 * monitored[183]
            + np.sqrt(monitored[183]) / (1 + np.sqrt(monitored[183]))
        )
        * monitored[162]
    )
    monitored[185] = np.exp(
        -(
            -0.3 * monitored[183]
            + np.sqrt(monitored[183]) / (1 + np.sqrt(monitored[183]))
        )
        * monitored[162]
    )
    monitored[186] = np.exp(
        -(
            -0.3 * monitored[183]
            + np.sqrt(monitored[183]) / (1 + np.sqrt(monitored[183]))
        )
        * monitored[162]
    )
    monitored[187] = (
        4
        * (-cao * monitored[164] + cai * np.exp(2 * monitored[273]) * monitored[184])
        * monitored[274]
        / (-1 + np.exp(2 * monitored[273]))
    )
    monitored[188] = (
        (-nao * monitored[166] + np.exp(monitored[273]) * monitored[185] * nai)
        * monitored[274]
        / (-1 + np.exp(monitored[273]))
    )
    monitored[189] = (
        (-ko * monitored[168] + np.exp(monitored[273]) * ki * monitored[186])
        * monitored[274]
        / (-1 + np.exp(monitored[273]))
    )
    monitored[190] = (1 - ICaL_fractionSS) * (
        (1 - monitored[178])
        * ((1 - nca_i) * monitored[146] + jca * monitored[152] * nca_i)
        * d
        * monitored[172]
        * monitored[187]
        + ((1 - nca_i) * monitored[155] + jca * monitored[157] * nca_i)
        * d
        * monitored[173]
        * monitored[178]
        * monitored[187]
    )
    monitored[191] = (1 - ICaL_fractionSS) * (
        (1 - monitored[178])
        * ((1 - nca_i) * monitored[146] + jca * monitored[152] * nca_i)
        * d
        * monitored[174]
        * monitored[188]
        + ((1 - nca_i) * monitored[155] + jca * monitored[157] * nca_i)
        * d
        * monitored[176]
        * monitored[178]
        * monitored[188]
    )
    monitored[192] = (1 - ICaL_fractionSS) * (
        (1 - monitored[178])
        * ((1 - nca_i) * monitored[146] + jca * monitored[152] * nca_i)
        * d
        * monitored[175]
        * monitored[189]
        + ((1 - nca_i) * monitored[155] + jca * monitored[157] * nca_i)
        * d
        * monitored[177]
        * monitored[178]
        * monitored[189]
    )
    monitored[193] = monitored[179] + monitored[190]
    monitored[194] = monitored[180] + monitored[191]
    monitored[195] = monitored[181] + monitored[192]

    # Expressions for the SERCA component
    monitored[196] = np.where(celltype == 1, 1.3, 1)
    monitored[197] = 0.005425 * cai * monitored[196] / (0.00092 + cai)
    monitored[198] = 0.01491875 * cai * monitored[196] / (0.00075 + cai)
    monitored[199] = 1.0 / (1 + KmCaMK / monitored[8])
    monitored[200] = 0.0003255 * cansr
    monitored[201] = Jup_b * (
        -monitored[200]
        + (1 - monitored[199]) * monitored[197]
        + monitored[198] * monitored[199]
    )

    # Expressions for the I_katp component
    monitored[202] = np.power(ko / K_o_n, 0.24)
    monitored[203] = 1.0 / (1 + (A_atp * A_atp) / (K_atp * K_atp))
    monitored[204] = (
        fkatp * gkatp * (-monitored[10] + v) * monitored[202] * monitored[203]
    )

    # Expressions for the INa component
    monitored[14] = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
    )
    monitored[15] = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    monitored[287] = (-m + monitored[14]) / monitored[15]
    monitored[16] = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
    )
    monitored[17] = np.where(
        v >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * v)
    )
    monitored[18] = np.where(
        v >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * v)),
        310000.0 * np.exp(0.3485 * v) + 2.7 * np.exp(0.079 * v),
    )
    monitored[19] = 1.0 / (monitored[17] + monitored[18])
    monitored[288] = (-h + monitored[16]) / monitored[19]
    monitored[20] = np.where(
        v >= -40,
        0,
        (37.78 + v)
        * (-25428.0 * np.exp(0.2444 * v) - 6.948e-06 * np.exp(-0.04391 * v))
        / (1 + 50262745825.95399 * np.exp(0.311 * v)),
    )
    monitored[21] = np.where(
        v >= -40,
        0.6 * np.exp(0.057 * v) / (1 + 0.040762203978366204 * np.exp(-0.1 * v)),
        0.02424
        * np.exp(-0.01052 * v)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * v)),
    )
    monitored[22] = monitored[16]
    monitored[23] = 1.0 / (monitored[20] + monitored[21])
    monitored[289] = (-j + monitored[22]) / monitored[23]
    monitored[24] = 1.0 / (
        (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
        * (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
    )
    monitored[290] = (-hp + monitored[24]) / monitored[19]
    monitored[25] = 1.46 * monitored[23]
    monitored[291] = (-jp + monitored[22]) / monitored[25]
    monitored[26] = 1.0 / (1 + KmCaMK / monitored[8])
    monitored[27] = (
        GNa
        * (m * m * m)
        * (-monitored[9] + v)
        * ((1 - monitored[26]) * h * j + hp * jp * monitored[26])
    )

    # Expressions for the L component
    monitored[127] = 1.0 / (1 + 0.000291579585635531 * np.exp(-0.18996960486322187 * v))
    monitored[128] = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    monitored[292] = (-mL + monitored[127]) / monitored[128]
    monitored[129] = 1.0 / (1 + 120578.15595522427 * np.exp(0.13354700854700854 * v))
    monitored[293] = (-hL + monitored[129]) / thL
    monitored[130] = 1.0 / (1 + 275969.2903869871 * np.exp(0.13354700854700854 * v))
    monitored[131] = 3 * thL
    monitored[294] = (-hLp + monitored[130]) / monitored[131]
    monitored[132] = np.where(celltype == 1, 0.6 * GNaL_b, GNaL_b)
    monitored[133] = 1.0 / (1 + KmCaMK / monitored[8])
    monitored[134] = (
        (-monitored[9] + v)
        * ((1 - monitored[133]) * hL + hLp * monitored[133])
        * mL
        * monitored[132]
    )

    # Expressions for the Ito component
    monitored[205] = 1.0 / (
        1
        + 2.6316508161673635
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    monitored[206] = 1.0515 / (
        1.0
        / (
            1.2089
            + 2.2621017070578837
            * np.exp(-0.03403513787634354 * EKshift - 0.03403513787634354 * v)
        )
        + 3.5
        / (
            1
            + 30.069572727397507
            * np.exp(0.03403513787634354 * EKshift + 0.03403513787634354 * v)
        )
    )
    monitored[295] = (-a + monitored[205]) / monitored[206]
    monitored[207] = 1.0 / (
        1
        + 2194.970764538301
        * np.exp(0.17510068289266328 * EKshift + 0.17510068289266328 * v)
    )
    monitored[208] = np.where(
        celltype == 1, 1 - 0.95 / (1 + np.exp(14 + EKshift / 5 + v / 5)), 1
    )
    monitored[209] = 4.562 + 1.0 / (
        0.3933 * np.exp(-1 - EKshift / 100 - v / 100)
        + 1.6300896349780942
        * np.exp(0.06027727546714889 * EKshift + 0.06027727546714889 * v)
    )
    monitored[210] = 23.62 + 1.0 / (
        0.00027617763953377436
        * np.exp(-0.01693480101608806 * EKshift - 0.01693480101608806 * v)
        + 0.024208962804604526
        * np.exp(0.12377769525931426 * EKshift + 0.12377769525931426 * v)
    )
    monitored[211] = monitored[208] * monitored[209]
    monitored[212] = monitored[208] * monitored[210]
    monitored[213] = 1.0 / (
        1
        + 0.24348537187522867
        * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
    )
    monitored[214] = 1 - monitored[213]
    monitored[296] = (-iF + monitored[207]) / monitored[211]
    monitored[297] = (-iS + monitored[207]) / monitored[212]
    monitored[215] = iF * monitored[213] + iS * monitored[214]
    monitored[216] = 1.0 / (
        1
        + 5.167428462230666
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    monitored[298] = (-ap + monitored[216]) / monitored[206]
    monitored[217] = 1.354 + 0.0001 / (
        2.6591269045230603e-05
        * np.exp(0.06293266205160478 * EKshift + 0.06293266205160478 * v)
        + 4.5541779737128264e24
        * np.exp(-4.642525533890436 * EKshift - 4.642525533890436 * v)
    )
    monitored[218] = 1 - 0.5 / (
        1 + 33.11545195869231 * np.exp(0.05 * EKshift + 0.05 * v)
    )
    monitored[219] = monitored[211] * monitored[217] * monitored[218]
    monitored[220] = monitored[212] * monitored[217] * monitored[218]
    monitored[299] = (-iFp + monitored[207]) / monitored[219]
    monitored[300] = (-iSp + monitored[207]) / monitored[220]
    monitored[221] = iFp * monitored[213] + iSp * monitored[214]
    monitored[222] = np.where(
        celltype == 1, 2 * Gto_b, np.where(celltype == 2, 2 * Gto_b, Gto_b)
    )
    monitored[223] = 1.0 / (1 + KmCaMK / monitored[8])
    monitored[224] = (
        (-monitored[10] + v)
        * (
            (1 - monitored[223]) * a * monitored[215]
            + ap * monitored[221] * monitored[223]
        )
        * monitored[222]
    )

    # Expressions for the IKr component
    monitored[225] = 0.1161 * np.exp(0.299 * monitored[273])
    monitored[226] = 0.2442 * np.exp(-1.604 * monitored[273])
    monitored[227] = 0.0578 * np.exp(0.971 * monitored[273])
    monitored[228] = 0.000349 * np.exp(-1.062 * monitored[273])
    monitored[229] = 0.2533 * np.exp(0.5953 * monitored[273])
    monitored[230] = 0.06525 * np.exp(-0.8209 * monitored[273])
    monitored[231] = 5.2e-05 * np.exp(1.525 * monitored[273])
    monitored[232] = (
        monitored[228]
        * monitored[230]
        * monitored[231]
        / (monitored[227] * monitored[229])
    )
    monitored[301] = C2 * monitored[226] - C3 * monitored[225]
    monitored[302] = beta_1 * C1 + C3 * monitored[225] - (alpha_1 + monitored[226]) * C2
    monitored[303] = (
        alpha_1 * C2
        + I * monitored[232]
        + O * monitored[228]
        - (beta_1 + monitored[227] + monitored[231]) * C1
    )
    monitored[304] = (
        C1 * monitored[227] + I * monitored[230] - (monitored[228] + monitored[229]) * O
    )
    monitored[305] = (
        C1 * monitored[231] + O * monitored[229] - (monitored[230] + monitored[232]) * I
    )
    monitored[233] = np.where(
        celltype == 1, 1.3 * GKr_b, np.where(celltype == 2, 0.8 * GKr_b, GKr_b)
    )
    monitored[234] = (
        np.sqrt(5) * np.sqrt(ko) * (-monitored[10] + v) * O * monitored[233] / 5
    )

    # Expressions for the IKs component
    monitored[235] = 1.0 / (1 + 0.27288596035656526 * np.exp(-0.11195700850873264 * v))
    monitored[236] = 817.3 + 1.0 / (
        0.003504067763074858 * np.exp(0.056179775280898875 * v)
        + 0.001292 * np.exp(-21 / 23 - v / 230)
    )
    monitored[306] = (-xs1 + monitored[235]) / monitored[236]
    monitored[237] = monitored[235]
    monitored[238] = 1.0 / (
        0.0022561357010639103 * np.exp(-v / 31) + 0.01 * np.exp(-5 / 2 + v / 20)
    )
    monitored[307] = (-xs2 + monitored[237]) / monitored[238]
    monitored[239] = 1 + 0.6 / (1 + 6.481821026062645e-07 * np.power(1.0 / cai, 1.4))
    monitored[240] = np.where(celltype == 1, 1.4 * GKs_b, GKs_b)
    monitored[241] = (-monitored[11] + v) * monitored[239] * monitored[240] * xs1 * xs2

    # Expressions for the IK1 component
    monitored[242] = 4.094 / (
        1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * monitored[10])
    )
    monitored[243] = (
        12.621629724407278 * np.exp(0.0674 * v - 0.0674 * monitored[10])
        + 1.1196358381249121e-16 * np.exp(0.0618 * v - 0.0618 * monitored[10])
    ) / (1 + 0.09883333819716558 * np.exp(0.1629 * monitored[10] - 0.1629 * v))
    monitored[244] = monitored[242] / (monitored[242] + monitored[243])
    monitored[245] = np.where(
        celltype == 1, 1.2 * GK1_b, np.where(celltype == 2, 1.3 * GK1_b, GK1_b)
    )
    monitored[246] = (
        np.sqrt(5)
        * np.sqrt(ko)
        * (-monitored[10] + v)
        * monitored[244]
        * monitored[245]
        / 5
    )

    # Expressions for the IKb component
    monitored[247] = 1.0 / (1 + 1.57503502085457 * np.exp(-0.04168907454423419 * v))
    monitored[248] = np.where(celltype == 1, 0.6 * GKb_b, GKb_b)
    monitored[249] = (-monitored[10] + v) * monitored[247] * monitored[248]

    # Expressions for the ICl component
    monitored[250] = Fjunc * GClCa * (-monitored[13] + v) / (1 + KdClCa / cass)
    monitored[251] = GClCa * (1 - Fjunc) * (-monitored[12] + v) / (1 + KdClCa / cai)
    monitored[252] = monitored[250] + monitored[251]
    monitored[253] = GClb * (-monitored[12] + v)

    # Expressions for the ICab component
    monitored[254] = (
        4
        * PCab
        * (-cao * monitored[164] + cai * np.exp(2 * monitored[273]) * monitored[184])
        * monitored[274]
        / (-1 + np.exp(2 * monitored[273]))
    )

    # Expressions for the Ryr component
    monitored[255] = 0.5 * bt
    monitored[256] = (
        -monitored[179]
        * monitored[255]
        / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    )
    monitored[257] = np.where(celltype == 2, 1.7 * monitored[256], monitored[256])
    monitored[258] = bt / (1 + 0.0123 / cajsr)
    monitored[259] = np.where(monitored[258] < 0.001, 0.001, monitored[258])
    monitored[308] = (-Jrel_np + monitored[257]) / monitored[259]
    monitored[260] = 1.25 * bt
    monitored[261] = 0.5 * monitored[260]
    monitored[262] = (
        -monitored[179]
        * monitored[261]
        / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    )
    monitored[263] = np.where(celltype == 2, 1.7 * monitored[262], monitored[262])
    monitored[264] = monitored[260] / (1 + 0.0123 / cajsr)
    monitored[265] = np.where(monitored[264] < 0.001, 0.001, monitored[264])
    monitored[309] = (-Jrel_p + monitored[263]) / monitored[265]
    monitored[266] = 1.0 / (1 + KmCaMK / monitored[8])
    monitored[267] = Jrel_b * ((1 - monitored[266]) * Jrel_np + Jrel_p * monitored[266])

    # Expressions for the Membrane component
    monitored[268] = np.where(
        np.logical_and(
            t >= i_Stim_Start,
            t
            - i_Stim_Start
            - i_Stim_Period * np.floor((t - i_Stim_Start) / i_Stim_Period)
            <= i_Stim_PulseDuration,
        ),
        i_Stim_Amplitude,
        0,
    )
    monitored[310] = (
        -monitored[102]
        - monitored[125]
        - monitored[126]
        - monitored[134]
        - monitored[135]
        - monitored[193]
        - monitored[194]
        - monitored[195]
        - monitored[204]
        - monitored[224]
        - monitored[234]
        - monitored[241]
        - monitored[246]
        - monitored[249]
        - monitored[252]
        - monitored[253]
        - monitored[254]
        - monitored[268]
        - monitored[27]
        - monitored[66]
    )

    # Expressions for the Intracellular ions component
    monitored[269] = np.where(celltype == 1, 1.3 * cmdnmax_b, cmdnmax_b)
    monitored[311] = monitored[136] * monitored[6] / monitored[3] + (
        -monitored[126]
        - monitored[134]
        - monitored[191]
        - monitored[27]
        - 3 * monitored[125]
        - 3 * monitored[66]
    ) * monitored[2] / (F * monitored[3])
    monitored[312] = -monitored[136] + (
        -monitored[180] - 3 * monitored[102]
    ) * monitored[2] / (F * monitored[6])
    monitored[313] = monitored[137] * monitored[6] / monitored[3] + (
        -monitored[192]
        - monitored[204]
        - monitored[224]
        - monitored[234]
        - monitored[241]
        - monitored[246]
        - monitored[249]
        - monitored[268]
        + 2 * monitored[125]
    ) * monitored[2] / (F * monitored[3])
    monitored[314] = -monitored[137] - monitored[181] * monitored[2] / (
        F * monitored[6]
    )
    monitored[315] = monitored[139] * monitored[6] / monitored[3] + (
        monitored[251] + monitored[253]
    ) * monitored[2] / (F * monitored[3])
    monitored[316] = -monitored[139] + monitored[250] * monitored[2] / (
        F * monitored[6]
    )
    monitored[270] = 1.0 / (
        1
        + kmcmdn * monitored[269] / ((kmcmdn + cai) * (kmcmdn + cai))
        + kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai))
    )
    monitored[317] = (
        monitored[138] * monitored[6] / monitored[3]
        - monitored[201] * monitored[4] / monitored[3]
        + (-monitored[135] - monitored[190] - monitored[254] + 2 * monitored[66])
        * monitored[2]
        / (2 * F * monitored[3])
    ) * monitored[270]
    monitored[271] = 1.0 / (
        1
        + BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass))
        + BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass))
    )
    monitored[318] = (
        -monitored[138]
        + monitored[267] * monitored[5] / monitored[6]
        + (-monitored[179] + 2 * monitored[102]) * monitored[2] / (2 * F * monitored[6])
    ) * monitored[271]
    monitored[319] = -monitored[275] * monitored[5] / monitored[4] + monitored[201]
    monitored[272] = 1.0 / (
        1 + csqnmax * kmcsqn / ((kmcsqn + cajsr) * (kmcsqn + cajsr))
    )
    monitored[320] = (-monitored[267] + monitored[275]) * monitored[272]

    # Return results
    return monitored


def forward_explicit_euler(states, t, dt, parameters):
    """
    Compute a forward step using the explicit Euler scheme to the\
        ToRORd_dyn_chloride_epi ODE
    """

    # Assign states
    assert len(states) == 45
    (
        CaMKt,
        d,
        ff,
        fs,
        fcaf,
        fcas,
        jca,
        ffp,
        fcafp,
        nca_ss,
        nca_i,
        m,
        h,
        j,
        hp,
        jp,
        mL,
        hL,
        hLp,
        a,
        iF,
        iS,
        ap,
        iFp,
        iSp,
        C3,
        C2,
        C1,
        O,
        I,
        xs1,
        xs2,
        Jrel_np,
        Jrel_p,
        v,
        nai,
        nass,
        ki,
        kss,
        cli,
        clss,
        cai,
        cass,
        cansr,
        cajsr,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    celltype = parameters[0]
    cao = parameters[1]
    clo = parameters[2]
    ko = parameters[3]
    nao = parameters[4]
    F = parameters[5]
    R = parameters[6]
    T = parameters[7]
    zca = parameters[8]
    zcl = parameters[9]
    zk = parameters[10]
    zna = parameters[11]
    L = parameters[12]
    rad = parameters[13]
    CaMKo = parameters[14]
    KmCaM = parameters[15]
    KmCaMK = parameters[16]
    aCaMK = parameters[17]
    bCaMK = parameters[18]
    PKNa = parameters[19]
    GNa = parameters[20]
    Gncx_b = parameters[21]
    INaCa_fractionSS = parameters[22]
    KmCaAct = parameters[23]
    kasymm = parameters[24]
    kcaoff = parameters[25]
    kcaon = parameters[26]
    kna1 = parameters[27]
    kna2 = parameters[28]
    kna3 = parameters[29]
    qca = parameters[30]
    qna = parameters[31]
    wca = parameters[32]
    wna = parameters[33]
    wnaca = parameters[34]
    H = parameters[35]
    Khp = parameters[36]
    Kki = parameters[37]
    Kko = parameters[38]
    Kmgatp = parameters[39]
    Knai0 = parameters[40]
    Knao0 = parameters[41]
    Knap = parameters[42]
    Kxkur = parameters[43]
    MgADP = parameters[44]
    MgATP = parameters[45]
    Pnak_b = parameters[46]
    delta = parameters[47]
    eP = parameters[48]
    k1m = parameters[49]
    k1p = parameters[50]
    k2m = parameters[51]
    k2p = parameters[52]
    k3m = parameters[53]
    k3p = parameters[54]
    k4m = parameters[55]
    k4p = parameters[56]
    PNab = parameters[57]
    GNaL_b = parameters[58]
    thL = parameters[59]
    GpCa = parameters[60]
    KmCap = parameters[61]
    tauCa = parameters[62]
    tauK = parameters[64]
    tauNa = parameters[65]
    Aff = parameters[66]
    ICaL_fractionSS = parameters[67]
    Kmn = parameters[68]
    PCa_b = parameters[69]
    dielConstant = parameters[70]
    k2n = parameters[71]
    offset = parameters[72]
    tjca = parameters[73]
    vShift = parameters[74]
    Jup_b = parameters[75]
    A_atp = parameters[76]
    K_atp = parameters[77]
    K_o_n = parameters[78]
    fkatp = parameters[79]
    gkatp = parameters[80]
    EKshift = parameters[81]
    Gto_b = parameters[82]
    GKr_b = parameters[83]
    alpha_1 = parameters[84]
    beta_1 = parameters[85]
    GKs_b = parameters[86]
    GK1_b = parameters[87]
    GKb_b = parameters[88]
    Fjunc = parameters[89]
    GClCa = parameters[90]
    GClb = parameters[91]
    KdClCa = parameters[92]
    PCab = parameters[93]
    Jrel_b = parameters[94]
    bt = parameters[95]
    cajsr_half = parameters[96]
    i_Stim_Amplitude = parameters[97]
    i_Stim_Period = parameters[99]
    i_Stim_PulseDuration = parameters[100]
    i_Stim_Start = parameters[101]
    BSLmax = parameters[102]
    BSRmax = parameters[103]
    KmBSL = parameters[104]
    KmBSR = parameters[105]
    cmdnmax_b = parameters[106]
    csqnmax = parameters[107]
    kmcmdn = parameters[108]
    kmcsqn = parameters[109]
    kmtrpn = parameters[110]
    trpnmax = parameters[111]

    # Expressions for the ToRORd dyn chloride component
    vfrt = F * v / (R * T)
    vffrt = (F * F) * v / (R * T)

    # Expressions for the Cell geometry component
    vcell = 3140.0 * L * (rad * rad)
    Ageo = 6.28 * (rad * rad) + 6.28 * L * rad
    Acap = 2 * Ageo
    vmyo = 0.68 * vcell
    vnsr = 0.0552 * vcell
    vjsr = 0.0048 * vcell
    vss = 0.02 * vcell

    # Expressions for the CaMK component
    CaMKb = CaMKo * (1 - CaMKt) / (1 + KmCaM / cass)
    CaMKa = CaMKb + CaMKt
    dCaMKt_dt = -bCaMK * CaMKt + aCaMK * (CaMKb + CaMKt) * CaMKb
    states[0] = dt * dCaMKt_dt + CaMKt

    # Expressions for the Reversal potentials component
    ENa = R * T * np.log(nao / nai) / (F * zna)
    EK = R * T * np.log(ko / ki) / (F * zk)
    EKs = R * T * np.log((ko + PKNa * nao) / (PKNa * nai + ki)) / (F * zk)
    ECl = R * T * np.log(clo / cli) / (F * zcl)
    EClss = R * T * np.log(clo / clss) / (F * zcl)

    # Expressions for the Ca component
    hca = np.exp(qca * vfrt)
    hna = np.exp(qna * vfrt)
    h1_i = 1 + (1 + hna) * nai / kna3
    h2_i = hna * nai / (kna3 * h1_i)
    h3_i = 1.0 / h1_i
    h4_i = 1 + (1 + nai / kna2) * nai / kna1
    h5_i = (nai * nai) / (kna1 * kna2 * h4_i)
    h6_i = 1.0 / h4_i
    h7_i = 1 + nao * (1 + 1.0 / hna) / kna3
    h8_i = nao / (kna3 * h7_i * hna)
    h9_i = 1.0 / h7_i
    h10_i = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    h11_i = (nao * nao) / (kna1 * kna2 * h10_i)
    h12_i = 1.0 / h10_i
    k1_i = cao * kcaon * h12_i
    k2_i = kcaoff
    k3p_i = wca * h9_i
    k3pp_i = wnaca * h8_i
    k3_i = k3p_i + k3pp_i
    k4p_i = wca * h3_i / hca
    k4pp_i = wnaca * h2_i
    k4_i = k4p_i + k4pp_i
    k5_i = kcaoff
    k6_i = kcaon * cai * h6_i
    k7_i = wna * h2_i * h5_i
    k8_i = wna * h11_i * h8_i
    x1_i = (k2_i + k3_i) * k5_i * k7_i + (k6_i + k7_i) * k2_i * k4_i
    x2_i = (k1_i + k8_i) * k4_i * k6_i + (k4_i + k5_i) * k1_i * k7_i
    x3_i = (k2_i + k3_i) * k6_i * k8_i + (k6_i + k7_i) * k1_i * k3_i
    x4_i = (k1_i + k8_i) * k3_i * k5_i + (k4_i + k5_i) * k2_i * k8_i
    E1_i = x1_i / (x1_i + x2_i + x3_i + x4_i)
    E2_i = x2_i / (x1_i + x2_i + x3_i + x4_i)
    E3_i = x3_i / (x1_i + x2_i + x3_i + x4_i)
    E4_i = x4_i / (x1_i + x2_i + x3_i + x4_i)
    allo_i = 1.0 / (1 + (KmCaAct * KmCaAct) / (cai * cai))
    JncxNa_i = E3_i * k4pp_i - E2_i * k3pp_i - 3 * E1_i * k8_i + 3 * E4_i * k7_i
    JncxCa_i = E2_i * k2_i - E1_i * k1_i
    Gncx = np.where(
        celltype == 1, 1.1 * Gncx_b, np.where(celltype == 2, 1.4 * Gncx_b, Gncx_b)
    )
    INaCa_i = (1 - INaCa_fractionSS) * (zca * JncxCa_i + zna * JncxNa_i) * Gncx * allo_i
    h1_ss = 1 + (1 + hna) * nass / kna3
    h2_ss = hna * nass / (kna3 * h1_ss)
    h3_ss = 1.0 / h1_ss
    h4_ss = 1 + (1 + nass / kna2) * nass / kna1
    h5_ss = (nass * nass) / (kna1 * kna2 * h4_ss)
    h6_ss = 1.0 / h4_ss
    h7_ss = 1 + nao * (1 + 1.0 / hna) / kna3
    h8_ss = nao / (kna3 * h7_ss * hna)
    h9_ss = 1.0 / h7_ss
    h10_ss = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    h11_ss = (nao * nao) / (kna1 * kna2 * h10_ss)
    h12_ss = 1.0 / h10_ss
    k1_ss = cao * kcaon * h12_ss
    k2_ss = kcaoff
    k3p_ss = wca * h9_ss
    k3pp_ss = wnaca * h8_ss
    k3_ss = k3p_ss + k3pp_ss
    k4p_ss = wca * h3_ss / hca
    k4pp_ss = wnaca * h2_ss
    k4_ss = k4p_ss + k4pp_ss
    k5_ss = kcaoff
    k6_ss = kcaon * cass * h6_ss
    k7_ss = wna * h2_ss * h5_ss
    k8_ss = wna * h11_ss * h8_ss
    x1_ss = (k2_ss + k3_ss) * k5_ss * k7_ss + (k6_ss + k7_ss) * k2_ss * k4_ss
    x2_ss = (k1_ss + k8_ss) * k4_ss * k6_ss + (k4_ss + k5_ss) * k1_ss * k7_ss
    x3_ss = (k2_ss + k3_ss) * k6_ss * k8_ss + (k6_ss + k7_ss) * k1_ss * k3_ss
    x4_ss = (k1_ss + k8_ss) * k3_ss * k5_ss + (k4_ss + k5_ss) * k2_ss * k8_ss
    E1_ss = x1_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E2_ss = x2_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E3_ss = x3_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E4_ss = x4_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    allo_ss = 1.0 / (1 + (KmCaAct * KmCaAct) / (cass * cass))
    JncxNa_ss = (
        E3_ss * k4pp_ss - E2_ss * k3pp_ss - 3 * E1_ss * k8_ss + 3 * E4_ss * k7_ss
    )
    JncxCa_ss = E2_ss * k2_ss - E1_ss * k1_ss
    INaCa_ss = INaCa_fractionSS * (zca * JncxCa_ss + zna * JncxNa_ss) * Gncx * allo_ss

    # Expressions for the K component
    Knai = Knai0 * np.exp(delta * vfrt / 3)
    Knao = Knao0 * np.exp((1 - delta) * vfrt / 3)
    P = eP / (1 + H / Khp + nai / Knap + ki / Kxkur)
    a1 = (
        k1p
        * (nai * nai * nai)
        / (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
            * (Knai * Knai * Knai)
        )
    )
    b1 = MgADP * k1m
    a2 = k2p
    b2 = (
        k2m
        * (nao * nao * nao)
        / (
            (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
            * (Knao * Knao * Knao)
        )
    )
    a3 = (
        k3p
        * (ko * ko)
        / (
            (Kko * Kko)
            * (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
        )
    )
    b3 = H * k3m * P / (1 + MgATP / Kmgatp)
    a4 = MgATP * k4p / (Kmgatp * (1 + MgATP / Kmgatp))
    b4 = (
        k4m
        * (ki * ki)
        / (
            (Kki * Kki)
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
        )
    )
    x1 = a1 * a2 * a4 + a1 * a2 * b3 + a2 * b3 * b4 + b2 * b3 * b4
    x2 = a1 * a2 * a3 + a2 * a3 * b4 + a3 * b1 * b4 + b1 * b2 * b4
    x3 = a2 * a3 * a4 + a3 * a4 * b1 + a4 * b1 * b2 + b1 * b2 * b3
    x4 = a1 * a3 * a4 + a1 * a4 * b2 + a1 * b2 * b3 + b2 * b3 * b4
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    JnakNa = -3 * E2 * b3 + 3 * E1 * a3
    JnakK = -2 * E3 * a1 + 2 * E4 * b1
    Pnak = np.where(
        celltype == 1, 0.9 * Pnak_b, np.where(celltype == 2, 0.7 * Pnak_b, Pnak_b)
    )
    INaK = (zk * JnakK + zna * JnakNa) * Pnak

    # Expressions for the b component
    INab = PNab * (-nao + np.exp(vfrt) * nai) * vffrt / (-1 + np.exp(vfrt))

    # Expressions for the IpCa component
    IpCa = GpCa * cai / (KmCap + cai)

    # Expressions for the Diff component
    JdiffNa = (-nai + nass) / tauNa
    JdiffK = (-ki + kss) / tauK
    Jdiff = (-cai + cass) / tauCa
    JdiffCl = (-cli + clss) / tauNa

    # Expressions for the Trans flux component
    Jtr = -cajsr / 60 + cansr / 60

    # Expressions for the ICaL component
    dss = np.where(v >= 31.4978, 1, 1.0763 * np.exp(-1.007 * np.exp(-0.0829 * v)))
    td = (
        0.6
        + offset
        + 1.0
        / (
            3.5254214873653824 * np.exp(0.09 * vShift + 0.09 * v)
            + 0.7408182206817179 * np.exp(-0.05 * vShift - 0.05 * v)
        )
    )
    dd_dt = (-d + dss) / td
    states[1] = dt * dd_dt + d
    fss = 1.0 / (1 + 199.86038496778565 * np.exp(0.27056277056277056 * v))
    tff = 7 + 1.0 / (0.0045 * np.exp(-2 - v / 10) + 0.0045 * np.exp(2 + v / 10))
    tfs = 1000 + 1.0 / (
        3.5e-05 * np.exp(-5 / 4 - v / 4) + 3.5e-05 * np.exp(5 / 6 + v / 6)
    )
    Afs = 1 - Aff
    dff_dt = (-ff + fss) / tff
    states[2] = dt * dff_dt + ff
    dfs_dt = (-fs + fss) / tfs
    states[3] = dt * dfs_dt + fs
    f = Aff * ff + Afs * fs
    fcass = fss
    tfcaf = 7 + 1.0 / (0.04 * np.exp(-4 / 7 + v / 7) + 0.04 * np.exp(4 / 7 - v / 7))
    tfcas = 100 + 1.0 / (0.00012 * np.exp(-v / 3) + 0.00012 * np.exp(v / 7))
    Afcaf = 0.3 + 0.6 / (1 + np.exp(-1 + v / 10))
    Afcas = 1 - Afcaf
    dfcaf_dt = (-fcaf + fcass) / tfcaf
    states[4] = dt * dfcaf_dt + fcaf
    dfcas_dt = (-fcas + fcass) / tfcas
    states[5] = dt * dfcas_dt + fcas
    fca = Afcaf * fcaf + Afcas * fcas
    jcass = 1.0 / (1.0 + 649.7401897235336 * np.exp(0.35821750967187277 * v))
    djca_dt = (-jca + jcass) / tjca
    states[6] = dt * djca_dt + jca
    tffp = 2.5 * tff
    dffp_dt = (-ffp + fss) / tffp
    states[7] = dt * dffp_dt + ffp
    fp = Aff * ffp + Afs * fs
    tfcafp = 2.5 * tfcaf
    dfcafp_dt = (-fcafp + fcass) / tfcafp
    states[8] = dt * dfcafp_dt + fcafp
    fcap = Afcaf * fcafp + Afcas * fcas
    km2n = jca
    anca_ss = 1.0 / (np.power(1 + Kmn / cass, 4) + k2n / km2n)
    dnca_ss_dt = k2n * anca_ss - km2n * nca_ss
    states[9] = dt * dnca_ss_dt + nca_ss
    Io = 0.0005 * clo + 0.0005 * ko + 0.0005 * nao + 0.002 * cao
    Iss = 0.0005 * clss + 0.0005 * kss + 0.0005 * nass + 0.002 * cass
    constA = 1820000.0 * np.power(T * dielConstant, -1.5)
    gamma_cass = np.exp(-4 * (-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_cao = np.exp(-4 * (-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    gamma_nass = np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_nao = np.exp(-(-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    gamma_kss = np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_ko = np.exp(-(-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    PhiCaL_ss = (
        4
        * (-cao * gamma_cao + cass * np.exp(2 * vfrt) * gamma_cass)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    PhiCaNa_ss = (
        (-nao * gamma_nao + np.exp(vfrt) * gamma_nass * nass)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    PhiCaK_ss = (
        (-ko * gamma_ko + np.exp(vfrt) * gamma_kss * kss) * vffrt / (-1 + np.exp(vfrt))
    )
    PCa = np.where(
        celltype == 1, 1.2 * PCa_b, np.where(celltype == 2, 2 * PCa_b, PCa_b)
    )
    PCap = 1.1 * PCa
    PCaNa = 0.00125 * PCa
    PCaK = 0.0003574 * PCa
    PCaNap = 0.00125 * PCap
    PCaKp = 0.0003574 * PCap
    fICaLp = 1.0 / (1 + KmCaMK / CaMKa)
    ICaL_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCa * PhiCaL_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCap * PhiCaL_ss * d * fICaLp
    )
    ICaNa_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaNa * PhiCaNa_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaNap * PhiCaNa_ss * d * fICaLp
    )
    ICaK_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaK * PhiCaK_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaKp * PhiCaK_ss * d * fICaLp
    )
    anca_i = 1.0 / (np.power(1 + Kmn / cai, 4) + k2n / km2n)
    dnca_i_dt = k2n * anca_i - km2n * nca_i
    states[10] = dt * dnca_i_dt + nca_i
    Ii = 0.0005 * cli + 0.0005 * ki + 0.0005 * nai + 0.002 * cai
    gamma_cai = np.exp(-4 * (-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    gamma_nai = np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    gamma_ki = np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    PhiCaL_i = (
        4
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    PhiCaNa_i = (
        (-nao * gamma_nao + np.exp(vfrt) * gamma_nai * nai)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    PhiCaK_i = (
        (-ko * gamma_ko + np.exp(vfrt) * gamma_ki * ki) * vffrt / (-1 + np.exp(vfrt))
    )
    ICaL_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCa * PhiCaL_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCap * PhiCaL_i * d * fICaLp
    )
    ICaNa_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaNa * PhiCaNa_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaNap * PhiCaNa_i * d * fICaLp
    )
    ICaK_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaK * PhiCaK_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaKp * PhiCaK_i * d * fICaLp
    )
    ICaL = ICaL_i + ICaL_ss
    ICaNa = ICaNa_i + ICaNa_ss
    ICaK = ICaK_i + ICaK_ss

    # Expressions for the SERCA component
    upScale = np.where(celltype == 1, 1.3, 1)
    Jupnp = 0.005425 * cai * upScale / (0.00092 + cai)
    Jupp = 0.01491875 * cai * upScale / (0.00075 + cai)
    fJupp = 1.0 / (1 + KmCaMK / CaMKa)
    Jleak = 0.0003255 * cansr
    Jup = Jup_b * (-Jleak + (1 - fJupp) * Jupnp + Jupp * fJupp)

    # Expressions for the I_katp component
    akik = np.power(ko / K_o_n, 0.24)
    bkik = 1.0 / (1 + (A_atp * A_atp) / (K_atp * K_atp))
    I_katp = fkatp * gkatp * (-EK + v) * akik * bkik

    # Expressions for the INa component
    mss = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
    )
    tm = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    dm_dt = (-m + mss) / tm
    states[11] = dt * dm_dt + m
    hss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
    )
    ah = np.where(
        v >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * v)
    )
    bh = np.where(
        v >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * v)),
        310000.0 * np.exp(0.3485 * v) + 2.7 * np.exp(0.079 * v),
    )
    th = 1.0 / (ah + bh)
    dh_dt = (-h + hss) / th
    states[12] = dt * dh_dt + h
    aj = np.where(
        v >= -40,
        0,
        (37.78 + v)
        * (-25428.0 * np.exp(0.2444 * v) - 6.948e-06 * np.exp(-0.04391 * v))
        / (1 + 50262745825.95399 * np.exp(0.311 * v)),
    )
    bj = np.where(
        v >= -40,
        0.6 * np.exp(0.057 * v) / (1 + 0.040762203978366204 * np.exp(-0.1 * v)),
        0.02424
        * np.exp(-0.01052 * v)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * v)),
    )
    jss = hss
    tj = 1.0 / (aj + bj)
    dj_dt = (-j + jss) / tj
    states[13] = dt * dj_dt + j
    hssp = 1.0 / (
        (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
        * (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
    )
    dhp_dt = (-hp + hssp) / th
    states[14] = dt * dhp_dt + hp
    tjp = 1.46 * tj
    djp_dt = (-jp + jss) / tjp
    states[15] = dt * djp_dt + jp
    fINap = 1.0 / (1 + KmCaMK / CaMKa)
    INa = GNa * (m * m * m) * (-ENa + v) * ((1 - fINap) * h * j + fINap * hp * jp)

    # Expressions for the L component
    mLss = 1.0 / (1 + 0.000291579585635531 * np.exp(-0.18996960486322187 * v))
    tmL = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    dmL_dt = (-mL + mLss) / tmL
    states[16] = dt * dmL_dt + mL
    hLss = 1.0 / (1 + 120578.15595522427 * np.exp(0.13354700854700854 * v))
    dhL_dt = (-hL + hLss) / thL
    states[17] = dt * dhL_dt + hL
    hLssp = 1.0 / (1 + 275969.2903869871 * np.exp(0.13354700854700854 * v))
    thLp = 3 * thL
    dhLp_dt = (-hLp + hLssp) / thLp
    states[18] = dt * dhLp_dt + hLp
    GNaL = np.where(celltype == 1, 0.6 * GNaL_b, GNaL_b)
    fINaLp = 1.0 / (1 + KmCaMK / CaMKa)
    INaL = (-ENa + v) * ((1 - fINaLp) * hL + fINaLp * hLp) * GNaL * mL

    # Expressions for the Ito component
    ass = 1.0 / (
        1
        + 2.6316508161673635
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    ta = 1.0515 / (
        1.0
        / (
            1.2089
            + 2.2621017070578837
            * np.exp(-0.03403513787634354 * EKshift - 0.03403513787634354 * v)
        )
        + 3.5
        / (
            1
            + 30.069572727397507
            * np.exp(0.03403513787634354 * EKshift + 0.03403513787634354 * v)
        )
    )
    da_dt = (-a + ass) / ta
    states[19] = dt * da_dt + a
    iss = 1.0 / (
        1
        + 2194.970764538301
        * np.exp(0.17510068289266328 * EKshift + 0.17510068289266328 * v)
    )
    delta_epi = np.where(
        celltype == 1, 1 - 0.95 / (1 + np.exp(14 + EKshift / 5 + v / 5)), 1
    )
    tiF_b = 4.562 + 1.0 / (
        0.3933 * np.exp(-1 - EKshift / 100 - v / 100)
        + 1.6300896349780942
        * np.exp(0.06027727546714889 * EKshift + 0.06027727546714889 * v)
    )
    tiS_b = 23.62 + 1.0 / (
        0.00027617763953377436
        * np.exp(-0.01693480101608806 * EKshift - 0.01693480101608806 * v)
        + 0.024208962804604526
        * np.exp(0.12377769525931426 * EKshift + 0.12377769525931426 * v)
    )
    tiF = delta_epi * tiF_b
    tiS = delta_epi * tiS_b
    AiF = 1.0 / (
        1
        + 0.24348537187522867
        * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
    )
    AiS = 1 - AiF
    diF_dt = (-iF + iss) / tiF
    states[20] = dt * diF_dt + iF
    diS_dt = (-iS + iss) / tiS
    states[21] = dt * diS_dt + iS
    i = AiF * iF + AiS * iS
    assp = 1.0 / (
        1
        + 5.167428462230666
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    dap_dt = (-ap + assp) / ta
    states[22] = dt * dap_dt + ap
    dti_develop = 1.354 + 0.0001 / (
        2.6591269045230603e-05
        * np.exp(0.06293266205160478 * EKshift + 0.06293266205160478 * v)
        + 4.5541779737128264e24
        * np.exp(-4.642525533890436 * EKshift - 4.642525533890436 * v)
    )
    dti_recover = 1 - 0.5 / (1 + 33.11545195869231 * np.exp(0.05 * EKshift + 0.05 * v))
    tiFp = dti_develop * dti_recover * tiF
    tiSp = dti_develop * dti_recover * tiS
    diFp_dt = (-iFp + iss) / tiFp
    states[23] = dt * diFp_dt + iFp
    diSp_dt = (-iSp + iss) / tiSp
    states[24] = dt * diSp_dt + iSp
    ip = AiF * iFp + AiS * iSp
    Gto = np.where(celltype == 1, 2 * Gto_b, np.where(celltype == 2, 2 * Gto_b, Gto_b))
    fItop = 1.0 / (1 + KmCaMK / CaMKa)
    Ito = (-EK + v) * ((1 - fItop) * a * i + ap * fItop * ip) * Gto

    # Expressions for the IKr component
    alpha = 0.1161 * np.exp(0.299 * vfrt)
    beta = 0.2442 * np.exp(-1.604 * vfrt)
    alpha_2 = 0.0578 * np.exp(0.971 * vfrt)
    beta_2 = 0.000349 * np.exp(-1.062 * vfrt)
    alpha_i = 0.2533 * np.exp(0.5953 * vfrt)
    beta_i = 0.06525 * np.exp(-0.8209 * vfrt)
    alpha_C2ToI = 5.2e-05 * np.exp(1.525 * vfrt)
    beta_ItoC2 = alpha_C2ToI * beta_2 * beta_i / (alpha_2 * alpha_i)
    dC3_dt = C2 * beta - C3 * alpha
    states[25] = dt * dC3_dt + C3
    dC2_dt = beta_1 * C1 + C3 * alpha - (alpha_1 + beta) * C2
    states[26] = dt * dC2_dt + C2
    dC1_dt = (
        alpha_1 * C2
        + I * beta_ItoC2
        + O * beta_2
        - (beta_1 + alpha_2 + alpha_C2ToI) * C1
    )
    states[27] = dt * dC1_dt + C1
    dO_dt = C1 * alpha_2 + I * beta_i - (alpha_i + beta_2) * O
    states[28] = dt * dO_dt + O
    dI_dt = C1 * alpha_C2ToI + O * alpha_i - (beta_ItoC2 + beta_i) * I
    states[29] = dt * dI_dt + I
    GKr = np.where(
        celltype == 1, 1.3 * GKr_b, np.where(celltype == 2, 0.8 * GKr_b, GKr_b)
    )
    IKr = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GKr * O / 5

    # Expressions for the IKs component
    xs1ss = 1.0 / (1 + 0.27288596035656526 * np.exp(-0.11195700850873264 * v))
    txs1 = 817.3 + 1.0 / (
        0.003504067763074858 * np.exp(0.056179775280898875 * v)
        + 0.001292 * np.exp(-21 / 23 - v / 230)
    )
    dxs1_dt = (-xs1 + xs1ss) / txs1
    states[30] = dt * dxs1_dt + xs1
    xs2ss = xs1ss
    txs2 = 1.0 / (
        0.0022561357010639103 * np.exp(-v / 31) + 0.01 * np.exp(-5 / 2 + v / 20)
    )
    dxs2_dt = (-xs2 + xs2ss) / txs2
    states[31] = dt * dxs2_dt + xs2
    KsCa = 1 + 0.6 / (1 + 6.481821026062645e-07 * np.power(1.0 / cai, 1.4))
    GKs = np.where(celltype == 1, 1.4 * GKs_b, GKs_b)
    IKs = (-EKs + v) * GKs * KsCa * xs1 * xs2

    # Expressions for the IK1 component
    aK1 = 4.094 / (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
    bK1 = (
        12.621629724407278 * np.exp(0.0674 * v - 0.0674 * EK)
        + 1.1196358381249121e-16 * np.exp(0.0618 * v - 0.0618 * EK)
    ) / (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
    K1ss = aK1 / (aK1 + bK1)
    GK1 = np.where(
        celltype == 1, 1.2 * GK1_b, np.where(celltype == 2, 1.3 * GK1_b, GK1_b)
    )
    IK1 = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GK1 * K1ss / 5

    # Expressions for the IKb component
    xkb = 1.0 / (1 + 1.57503502085457 * np.exp(-0.04168907454423419 * v))
    GKb = np.where(celltype == 1, 0.6 * GKb_b, GKb_b)
    IKb = (-EK + v) * GKb * xkb

    # Expressions for the ICl component
    IClCa_junc = Fjunc * GClCa * (-EClss + v) / (1 + KdClCa / cass)
    IClCa_sl = GClCa * (1 - Fjunc) * (-ECl + v) / (1 + KdClCa / cai)
    IClCa = IClCa_junc + IClCa_sl
    IClb = GClb * (-ECl + v)

    # Expressions for the ICab component
    ICab = (
        4
        * PCab
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )

    # Expressions for the Ryr component
    a_rel = 0.5 * bt
    Jrel_inf_b = -ICaL_ss * a_rel / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    Jrel_inf = np.where(celltype == 2, 1.7 * Jrel_inf_b, Jrel_inf_b)
    tau_rel_b = bt / (1 + 0.0123 / cajsr)
    tau_rel = np.where(tau_rel_b < 0.001, 0.001, tau_rel_b)
    dJrel_np_dt = (-Jrel_np + Jrel_inf) / tau_rel
    states[32] = dt * dJrel_np_dt + Jrel_np
    btp = 1.25 * bt
    a_relp = 0.5 * btp
    Jrel_infp_b = -ICaL_ss * a_relp / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    Jrel_infp = np.where(celltype == 2, 1.7 * Jrel_infp_b, Jrel_infp_b)
    tau_relp_b = btp / (1 + 0.0123 / cajsr)
    tau_relp = np.where(tau_relp_b < 0.001, 0.001, tau_relp_b)
    dJrel_p_dt = (-Jrel_p + Jrel_infp) / tau_relp
    states[33] = dt * dJrel_p_dt + Jrel_p
    fJrelp = 1.0 / (1 + KmCaMK / CaMKa)
    Jrel = Jrel_b * ((1 - fJrelp) * Jrel_np + Jrel_p * fJrelp)

    # Expressions for the Membrane component
    Istim = np.where(
        np.logical_and(
            t >= i_Stim_Start,
            t
            - i_Stim_Start
            - i_Stim_Period * np.floor((t - i_Stim_Start) / i_Stim_Period)
            <= i_Stim_PulseDuration,
        ),
        i_Stim_Amplitude,
        0,
    )
    dv_dt = (
        -ICaK
        - ICaL
        - ICaNa
        - ICab
        - IClCa
        - IClb
        - IK1
        - IKb
        - IKr
        - IKs
        - INa
        - INaCa_i
        - INaCa_ss
        - INaK
        - INaL
        - INab
        - I_katp
        - IpCa
        - Istim
        - Ito
    )
    states[34] = dt * dv_dt + v

    # Expressions for the Intracellular ions component
    cmdnmax = np.where(celltype == 1, 1.3 * cmdnmax_b, cmdnmax_b)
    dnai_dt = JdiffNa * vss / vmyo + (
        -ICaNa_i - INa - INaL - INab - 3 * INaCa_i - 3 * INaK
    ) * Acap / (F * vmyo)
    states[35] = dt * dnai_dt + nai
    dnass_dt = -JdiffNa + (-ICaNa_ss - 3 * INaCa_ss) * Acap / (F * vss)
    states[36] = dt * dnass_dt + nass
    dki_dt = JdiffK * vss / vmyo + (
        -ICaK_i - IK1 - IKb - IKr - IKs - I_katp - Istim - Ito + 2 * INaK
    ) * Acap / (F * vmyo)
    states[37] = dt * dki_dt + ki
    dkss_dt = -JdiffK - Acap * ICaK_ss / (F * vss)
    states[38] = dt * dkss_dt + kss
    dcli_dt = JdiffCl * vss / vmyo + (IClCa_sl + IClb) * Acap / (F * vmyo)
    states[39] = dt * dcli_dt + cli
    dclss_dt = -JdiffCl + Acap * IClCa_junc / (F * vss)
    states[40] = dt * dclss_dt + clss
    Bcai = 1.0 / (
        1
        + kmcmdn * cmdnmax / ((kmcmdn + cai) * (kmcmdn + cai))
        + kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai))
    )
    dcai_dt = (
        Jdiff * vss / vmyo
        - Jup * vnsr / vmyo
        + (-ICaL_i - ICab - IpCa + 2 * INaCa_i) * Acap / (2 * F * vmyo)
    ) * Bcai
    states[41] = dt * dcai_dt + cai
    Bcass = 1.0 / (
        1
        + BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass))
        + BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass))
    )
    dcass_dt = (
        -Jdiff + Jrel * vjsr / vss + (-ICaL_ss + 2 * INaCa_ss) * Acap / (2 * F * vss)
    ) * Bcass
    states[42] = dt * dcass_dt + cass
    dcansr_dt = -Jtr * vjsr / vnsr + Jup
    states[43] = dt * dcansr_dt + cansr
    Bcajsr = 1.0 / (1 + csqnmax * kmcsqn / ((kmcsqn + cajsr) * (kmcsqn + cajsr)))
    dcajsr_dt = (-Jrel + Jtr) * Bcajsr
    states[44] = dt * dcajsr_dt + cajsr

    # Return results
    return states


def forward_generalized_rush_larsen(states, t, dt, parameters):
    """
    Compute a forward step using the generalised Rush-Larsen (GRL1) scheme to\
        the ToRORd_dyn_chloride_epi ODE
    """

    # Assign states
    assert len(states) == 45
    (
        CaMKt,
        d,
        ff,
        fs,
        fcaf,
        fcas,
        jca,
        ffp,
        fcafp,
        nca_ss,
        nca_i,
        m,
        h,
        j,
        hp,
        jp,
        mL,
        hL,
        hLp,
        a,
        iF,
        iS,
        ap,
        iFp,
        iSp,
        C3,
        C2,
        C1,
        O,
        I,
        xs1,
        xs2,
        Jrel_np,
        Jrel_p,
        v,
        nai,
        nass,
        ki,
        kss,
        cli,
        clss,
        cai,
        cass,
        cansr,
        cajsr,
    ) = states

    # Assign parameters
    assert len(parameters) == 112
    celltype = parameters[0]
    cao = parameters[1]
    clo = parameters[2]
    ko = parameters[3]
    nao = parameters[4]
    F = parameters[5]
    R = parameters[6]
    T = parameters[7]
    zca = parameters[8]
    zcl = parameters[9]
    zk = parameters[10]
    zna = parameters[11]
    L = parameters[12]
    rad = parameters[13]
    CaMKo = parameters[14]
    KmCaM = parameters[15]
    KmCaMK = parameters[16]
    aCaMK = parameters[17]
    bCaMK = parameters[18]
    PKNa = parameters[19]
    GNa = parameters[20]
    Gncx_b = parameters[21]
    INaCa_fractionSS = parameters[22]
    KmCaAct = parameters[23]
    kasymm = parameters[24]
    kcaoff = parameters[25]
    kcaon = parameters[26]
    kna1 = parameters[27]
    kna2 = parameters[28]
    kna3 = parameters[29]
    qca = parameters[30]
    qna = parameters[31]
    wca = parameters[32]
    wna = parameters[33]
    wnaca = parameters[34]
    H = parameters[35]
    Khp = parameters[36]
    Kki = parameters[37]
    Kko = parameters[38]
    Kmgatp = parameters[39]
    Knai0 = parameters[40]
    Knao0 = parameters[41]
    Knap = parameters[42]
    Kxkur = parameters[43]
    MgADP = parameters[44]
    MgATP = parameters[45]
    Pnak_b = parameters[46]
    delta = parameters[47]
    eP = parameters[48]
    k1m = parameters[49]
    k1p = parameters[50]
    k2m = parameters[51]
    k2p = parameters[52]
    k3m = parameters[53]
    k3p = parameters[54]
    k4m = parameters[55]
    k4p = parameters[56]
    PNab = parameters[57]
    GNaL_b = parameters[58]
    thL = parameters[59]
    GpCa = parameters[60]
    KmCap = parameters[61]
    tauCa = parameters[62]
    tauK = parameters[64]
    tauNa = parameters[65]
    Aff = parameters[66]
    ICaL_fractionSS = parameters[67]
    Kmn = parameters[68]
    PCa_b = parameters[69]
    dielConstant = parameters[70]
    k2n = parameters[71]
    offset = parameters[72]
    tjca = parameters[73]
    vShift = parameters[74]
    Jup_b = parameters[75]
    A_atp = parameters[76]
    K_atp = parameters[77]
    K_o_n = parameters[78]
    fkatp = parameters[79]
    gkatp = parameters[80]
    EKshift = parameters[81]
    Gto_b = parameters[82]
    GKr_b = parameters[83]
    alpha_1 = parameters[84]
    beta_1 = parameters[85]
    GKs_b = parameters[86]
    GK1_b = parameters[87]
    GKb_b = parameters[88]
    Fjunc = parameters[89]
    GClCa = parameters[90]
    GClb = parameters[91]
    KdClCa = parameters[92]
    PCab = parameters[93]
    Jrel_b = parameters[94]
    bt = parameters[95]
    cajsr_half = parameters[96]
    i_Stim_Amplitude = parameters[97]
    i_Stim_Period = parameters[99]
    i_Stim_PulseDuration = parameters[100]
    i_Stim_Start = parameters[101]
    BSLmax = parameters[102]
    BSRmax = parameters[103]
    KmBSL = parameters[104]
    KmBSR = parameters[105]
    cmdnmax_b = parameters[106]
    csqnmax = parameters[107]
    kmcmdn = parameters[108]
    kmcsqn = parameters[109]
    kmtrpn = parameters[110]
    trpnmax = parameters[111]

    # Expressions for the ToRORd dyn chloride component
    vfrt = F * v / (R * T)
    vffrt = (F * F) * v / (R * T)

    # Expressions for the Cell geometry component
    vcell = 3140.0 * L * (rad * rad)
    Ageo = 6.28 * (rad * rad) + 6.28 * L * rad
    Acap = 2 * Ageo
    vmyo = 0.68 * vcell
    vnsr = 0.0552 * vcell
    vjsr = 0.0048 * vcell
    vss = 0.02 * vcell

    # Expressions for the CaMK component
    CaMKb = CaMKo * (1 - CaMKt) / (1 + KmCaM / cass)
    CaMKa = CaMKb + CaMKt
    dCaMKt_dt = -bCaMK * CaMKt + aCaMK * (CaMKb + CaMKt) * CaMKb
    dCaMKb_dCaMKt = -CaMKo / (1 + KmCaM / cass)
    dCaMKt_dt_linearized = (
        -bCaMK
        + aCaMK * (1 + dCaMKb_dCaMKt) * CaMKb
        + aCaMK * (CaMKb + CaMKt) * dCaMKb_dCaMKt
    )
    states[0] = CaMKt + np.where(
        np.abs(dCaMKt_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dCaMKt_dt_linearized)) * dCaMKt_dt / dCaMKt_dt_linearized,
        dt * dCaMKt_dt,
    )

    # Expressions for the Reversal potentials component
    ENa = R * T * np.log(nao / nai) / (F * zna)
    EK = R * T * np.log(ko / ki) / (F * zk)
    EKs = R * T * np.log((ko + PKNa * nao) / (PKNa * nai + ki)) / (F * zk)
    ECl = R * T * np.log(clo / cli) / (F * zcl)
    EClss = R * T * np.log(clo / clss) / (F * zcl)

    # Expressions for the Ca component
    hca = np.exp(qca * vfrt)
    hna = np.exp(qna * vfrt)
    h1_i = 1 + (1 + hna) * nai / kna3
    h2_i = hna * nai / (kna3 * h1_i)
    h3_i = 1.0 / h1_i
    h4_i = 1 + (1 + nai / kna2) * nai / kna1
    h5_i = (nai * nai) / (kna1 * kna2 * h4_i)
    h6_i = 1.0 / h4_i
    h7_i = 1 + nao * (1 + 1.0 / hna) / kna3
    h8_i = nao / (kna3 * h7_i * hna)
    h9_i = 1.0 / h7_i
    h10_i = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    h11_i = (nao * nao) / (kna1 * kna2 * h10_i)
    h12_i = 1.0 / h10_i
    k1_i = cao * kcaon * h12_i
    k2_i = kcaoff
    k3p_i = wca * h9_i
    k3pp_i = wnaca * h8_i
    k3_i = k3p_i + k3pp_i
    k4p_i = wca * h3_i / hca
    k4pp_i = wnaca * h2_i
    k4_i = k4p_i + k4pp_i
    k5_i = kcaoff
    k6_i = kcaon * cai * h6_i
    k7_i = wna * h2_i * h5_i
    k8_i = wna * h11_i * h8_i
    x1_i = (k2_i + k3_i) * k5_i * k7_i + (k6_i + k7_i) * k2_i * k4_i
    x2_i = (k1_i + k8_i) * k4_i * k6_i + (k4_i + k5_i) * k1_i * k7_i
    x3_i = (k2_i + k3_i) * k6_i * k8_i + (k6_i + k7_i) * k1_i * k3_i
    x4_i = (k1_i + k8_i) * k3_i * k5_i + (k4_i + k5_i) * k2_i * k8_i
    E1_i = x1_i / (x1_i + x2_i + x3_i + x4_i)
    E2_i = x2_i / (x1_i + x2_i + x3_i + x4_i)
    E3_i = x3_i / (x1_i + x2_i + x3_i + x4_i)
    E4_i = x4_i / (x1_i + x2_i + x3_i + x4_i)
    allo_i = 1.0 / (1 + (KmCaAct * KmCaAct) / (cai * cai))
    JncxNa_i = E3_i * k4pp_i - E2_i * k3pp_i - 3 * E1_i * k8_i + 3 * E4_i * k7_i
    JncxCa_i = E2_i * k2_i - E1_i * k1_i
    Gncx = np.where(
        celltype == 1, 1.1 * Gncx_b, np.where(celltype == 2, 1.4 * Gncx_b, Gncx_b)
    )
    INaCa_i = (1 - INaCa_fractionSS) * (zca * JncxCa_i + zna * JncxNa_i) * Gncx * allo_i
    h1_ss = 1 + (1 + hna) * nass / kna3
    h2_ss = hna * nass / (kna3 * h1_ss)
    h3_ss = 1.0 / h1_ss
    h4_ss = 1 + (1 + nass / kna2) * nass / kna1
    h5_ss = (nass * nass) / (kna1 * kna2 * h4_ss)
    h6_ss = 1.0 / h4_ss
    h7_ss = 1 + nao * (1 + 1.0 / hna) / kna3
    h8_ss = nao / (kna3 * h7_ss * hna)
    h9_ss = 1.0 / h7_ss
    h10_ss = 1 + kasymm + nao * (1 + nao / kna2) / kna1
    h11_ss = (nao * nao) / (kna1 * kna2 * h10_ss)
    h12_ss = 1.0 / h10_ss
    k1_ss = cao * kcaon * h12_ss
    k2_ss = kcaoff
    k3p_ss = wca * h9_ss
    k3pp_ss = wnaca * h8_ss
    k3_ss = k3p_ss + k3pp_ss
    k4p_ss = wca * h3_ss / hca
    k4pp_ss = wnaca * h2_ss
    k4_ss = k4p_ss + k4pp_ss
    k5_ss = kcaoff
    k6_ss = kcaon * cass * h6_ss
    k7_ss = wna * h2_ss * h5_ss
    k8_ss = wna * h11_ss * h8_ss
    x1_ss = (k2_ss + k3_ss) * k5_ss * k7_ss + (k6_ss + k7_ss) * k2_ss * k4_ss
    x2_ss = (k1_ss + k8_ss) * k4_ss * k6_ss + (k4_ss + k5_ss) * k1_ss * k7_ss
    x3_ss = (k2_ss + k3_ss) * k6_ss * k8_ss + (k6_ss + k7_ss) * k1_ss * k3_ss
    x4_ss = (k1_ss + k8_ss) * k3_ss * k5_ss + (k4_ss + k5_ss) * k2_ss * k8_ss
    E1_ss = x1_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E2_ss = x2_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E3_ss = x3_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    E4_ss = x4_ss / (x1_ss + x2_ss + x3_ss + x4_ss)
    allo_ss = 1.0 / (1 + (KmCaAct * KmCaAct) / (cass * cass))
    JncxNa_ss = (
        E3_ss * k4pp_ss - E2_ss * k3pp_ss - 3 * E1_ss * k8_ss + 3 * E4_ss * k7_ss
    )
    JncxCa_ss = E2_ss * k2_ss - E1_ss * k1_ss
    INaCa_ss = INaCa_fractionSS * (zca * JncxCa_ss + zna * JncxNa_ss) * Gncx * allo_ss

    # Expressions for the K component
    Knai = Knai0 * np.exp(delta * vfrt / 3)
    Knao = Knao0 * np.exp((1 - delta) * vfrt / 3)
    P = eP / (1 + H / Khp + nai / Knap + ki / Kxkur)
    a1 = (
        k1p
        * (nai * nai * nai)
        / (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
            * (Knai * Knai * Knai)
        )
    )
    b1 = MgADP * k1m
    a2 = k2p
    b2 = (
        k2m
        * (nao * nao * nao)
        / (
            (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
            * (Knao * Knao * Knao)
        )
    )
    a3 = (
        k3p
        * (ko * ko)
        / (
            (Kko * Kko)
            * (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
        )
    )
    b3 = H * k3m * P / (1 + MgATP / Kmgatp)
    a4 = MgATP * k4p / (Kmgatp * (1 + MgATP / Kmgatp))
    b4 = (
        k4m
        * (ki * ki)
        / (
            (Kki * Kki)
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
        )
    )
    x1 = a1 * a2 * a4 + a1 * a2 * b3 + a2 * b3 * b4 + b2 * b3 * b4
    x2 = a1 * a2 * a3 + a2 * a3 * b4 + a3 * b1 * b4 + b1 * b2 * b4
    x3 = a2 * a3 * a4 + a3 * a4 * b1 + a4 * b1 * b2 + b1 * b2 * b3
    x4 = a1 * a3 * a4 + a1 * a4 * b2 + a1 * b2 * b3 + b2 * b3 * b4
    E1 = x1 / (x1 + x2 + x3 + x4)
    E2 = x2 / (x1 + x2 + x3 + x4)
    E3 = x3 / (x1 + x2 + x3 + x4)
    E4 = x4 / (x1 + x2 + x3 + x4)
    JnakNa = -3 * E2 * b3 + 3 * E1 * a3
    JnakK = -2 * E3 * a1 + 2 * E4 * b1
    Pnak = np.where(
        celltype == 1, 0.9 * Pnak_b, np.where(celltype == 2, 0.7 * Pnak_b, Pnak_b)
    )
    INaK = (zk * JnakK + zna * JnakNa) * Pnak

    # Expressions for the b component
    INab = PNab * (-nao + np.exp(vfrt) * nai) * vffrt / (-1 + np.exp(vfrt))

    # Expressions for the IpCa component
    IpCa = GpCa * cai / (KmCap + cai)

    # Expressions for the Diff component
    JdiffNa = (-nai + nass) / tauNa
    JdiffK = (-ki + kss) / tauK
    Jdiff = (-cai + cass) / tauCa
    JdiffCl = (-cli + clss) / tauNa

    # Expressions for the Trans flux component
    Jtr = -cajsr / 60 + cansr / 60

    # Expressions for the ICaL component
    dss = np.where(v >= 31.4978, 1, 1.0763 * np.exp(-1.007 * np.exp(-0.0829 * v)))
    td = (
        0.6
        + offset
        + 1.0
        / (
            3.5254214873653824 * np.exp(0.09 * vShift + 0.09 * v)
            + 0.7408182206817179 * np.exp(-0.05 * vShift - 0.05 * v)
        )
    )
    dd_dt = (-d + dss) / td
    dd_dt_linearized = -1 / td
    states[1] = (-1 + np.exp(dt * dd_dt_linearized)) * dd_dt / dd_dt_linearized + d
    fss = 1.0 / (1 + 199.86038496778565 * np.exp(0.27056277056277056 * v))
    tff = 7 + 1.0 / (0.0045 * np.exp(-2 - v / 10) + 0.0045 * np.exp(2 + v / 10))
    tfs = 1000 + 1.0 / (
        3.5e-05 * np.exp(-5 / 4 - v / 4) + 3.5e-05 * np.exp(5 / 6 + v / 6)
    )
    Afs = 1 - Aff
    dff_dt = (-ff + fss) / tff
    dff_dt_linearized = -1 / tff
    states[2] = (-1 + np.exp(dt * dff_dt_linearized)) * dff_dt / dff_dt_linearized + ff
    dfs_dt = (-fs + fss) / tfs
    dfs_dt_linearized = -1 / tfs
    states[3] = (-1 + np.exp(dt * dfs_dt_linearized)) * dfs_dt / dfs_dt_linearized + fs
    f = Aff * ff + Afs * fs
    fcass = fss
    tfcaf = 7 + 1.0 / (0.04 * np.exp(-4 / 7 + v / 7) + 0.04 * np.exp(4 / 7 - v / 7))
    tfcas = 100 + 1.0 / (0.00012 * np.exp(-v / 3) + 0.00012 * np.exp(v / 7))
    Afcaf = 0.3 + 0.6 / (1 + np.exp(-1 + v / 10))
    Afcas = 1 - Afcaf
    dfcaf_dt = (-fcaf + fcass) / tfcaf
    dfcaf_dt_linearized = -1 / tfcaf
    states[4] = (
        -1 + np.exp(dt * dfcaf_dt_linearized)
    ) * dfcaf_dt / dfcaf_dt_linearized + fcaf
    dfcas_dt = (-fcas + fcass) / tfcas
    dfcas_dt_linearized = -1 / tfcas
    states[5] = (
        -1 + np.exp(dt * dfcas_dt_linearized)
    ) * dfcas_dt / dfcas_dt_linearized + fcas
    fca = Afcaf * fcaf + Afcas * fcas
    jcass = 1.0 / (1.0 + 649.7401897235336 * np.exp(0.35821750967187277 * v))
    djca_dt = (-jca + jcass) / tjca
    djca_dt_linearized = -1 / tjca
    states[6] = (
        -1 + np.exp(dt * djca_dt_linearized)
    ) * djca_dt / djca_dt_linearized + jca
    tffp = 2.5 * tff
    dffp_dt = (-ffp + fss) / tffp
    dffp_dt_linearized = -1 / tffp
    states[7] = (
        -1 + np.exp(dt * dffp_dt_linearized)
    ) * dffp_dt / dffp_dt_linearized + ffp
    fp = Aff * ffp + Afs * fs
    tfcafp = 2.5 * tfcaf
    dfcafp_dt = (-fcafp + fcass) / tfcafp
    dfcafp_dt_linearized = -1 / tfcafp
    states[8] = (
        -1 + np.exp(dt * dfcafp_dt_linearized)
    ) * dfcafp_dt / dfcafp_dt_linearized + fcafp
    fcap = Afcaf * fcafp + Afcas * fcas
    km2n = jca
    anca_ss = 1.0 / (np.power(1 + Kmn / cass, 4) + k2n / km2n)
    dnca_ss_dt = k2n * anca_ss - km2n * nca_ss
    dnca_ss_dt_linearized = -km2n
    states[9] = (
        np.where(
            np.abs(dnca_ss_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dnca_ss_dt_linearized))
            * dnca_ss_dt
            / dnca_ss_dt_linearized,
            dt * dnca_ss_dt,
        )
        + nca_ss
    )
    Io = 0.0005 * clo + 0.0005 * ko + 0.0005 * nao + 0.002 * cao
    Iss = 0.0005 * clss + 0.0005 * kss + 0.0005 * nass + 0.002 * cass
    constA = 1820000.0 * np.power(T * dielConstant, -1.5)
    gamma_cass = np.exp(-4 * (-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_cao = np.exp(-4 * (-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    gamma_nass = np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_nao = np.exp(-(-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    gamma_kss = np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    gamma_ko = np.exp(-(-0.3 * Io + np.sqrt(Io) / (1 + np.sqrt(Io))) * constA)
    PhiCaL_ss = (
        4
        * (-cao * gamma_cao + cass * np.exp(2 * vfrt) * gamma_cass)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    PhiCaNa_ss = (
        (-nao * gamma_nao + np.exp(vfrt) * gamma_nass * nass)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    PhiCaK_ss = (
        (-ko * gamma_ko + np.exp(vfrt) * gamma_kss * kss) * vffrt / (-1 + np.exp(vfrt))
    )
    PCa = np.where(
        celltype == 1, 1.2 * PCa_b, np.where(celltype == 2, 2 * PCa_b, PCa_b)
    )
    PCap = 1.1 * PCa
    PCaNa = 0.00125 * PCa
    PCaK = 0.0003574 * PCa
    PCaNap = 0.00125 * PCap
    PCaKp = 0.0003574 * PCap
    fICaLp = 1.0 / (1 + KmCaMK / CaMKa)
    ICaL_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCa * PhiCaL_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCap * PhiCaL_ss * d * fICaLp
    )
    ICaNa_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaNa * PhiCaNa_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaNap * PhiCaNa_ss * d * fICaLp
    )
    ICaK_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaK * PhiCaK_ss * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaKp * PhiCaK_ss * d * fICaLp
    )
    anca_i = 1.0 / (np.power(1 + Kmn / cai, 4) + k2n / km2n)
    dnca_i_dt = k2n * anca_i - km2n * nca_i
    dnca_i_dt_linearized = -km2n
    states[10] = (
        np.where(
            np.abs(dnca_i_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dnca_i_dt_linearized)) * dnca_i_dt / dnca_i_dt_linearized,
            dt * dnca_i_dt,
        )
        + nca_i
    )
    Ii = 0.0005 * cli + 0.0005 * ki + 0.0005 * nai + 0.002 * cai
    gamma_cai = np.exp(-4 * (-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    gamma_nai = np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    gamma_ki = np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    PhiCaL_i = (
        4
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    PhiCaNa_i = (
        (-nao * gamma_nao + np.exp(vfrt) * gamma_nai * nai)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    PhiCaK_i = (
        (-ko * gamma_ko + np.exp(vfrt) * gamma_ki * ki) * vffrt / (-1 + np.exp(vfrt))
    )
    ICaL_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCa * PhiCaL_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCap * PhiCaL_i * d * fICaLp
    )
    ICaNa_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaNa * PhiCaNa_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaNap * PhiCaNa_i * d * fICaLp
    )
    ICaK_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaK * PhiCaK_i * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaKp * PhiCaK_i * d * fICaLp
    )
    ICaL = ICaL_i + ICaL_ss
    ICaNa = ICaNa_i + ICaNa_ss
    ICaK = ICaK_i + ICaK_ss

    # Expressions for the SERCA component
    upScale = np.where(celltype == 1, 1.3, 1)
    Jupnp = 0.005425 * cai * upScale / (0.00092 + cai)
    Jupp = 0.01491875 * cai * upScale / (0.00075 + cai)
    fJupp = 1.0 / (1 + KmCaMK / CaMKa)
    Jleak = 0.0003255 * cansr
    Jup = Jup_b * (-Jleak + (1 - fJupp) * Jupnp + Jupp * fJupp)

    # Expressions for the I_katp component
    akik = np.power(ko / K_o_n, 0.24)
    bkik = 1.0 / (1 + (A_atp * A_atp) / (K_atp * K_atp))
    I_katp = fkatp * gkatp * (-EK + v) * akik * bkik

    # Expressions for the INa component
    mss = 1.0 / (
        (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
        * (1 + 0.0018422115811651339 * np.exp(-0.1107419712070875 * v))
    )
    tm = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    dm_dt = (-m + mss) / tm
    dm_dt_linearized = -1 / tm
    states[11] = (-1 + np.exp(dt * dm_dt_linearized)) * dm_dt / dm_dt_linearized + m
    hss = 1.0 / (
        (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
        * (1 + 15212.593285654404 * np.exp(0.13458950201884254 * v))
    )
    ah = np.where(
        v >= -40, 0, 4.4312679295805147e-07 * np.exp(-0.14705882352941177 * v)
    )
    bh = np.where(
        v >= -40,
        0.77 / (0.13 + 0.049758141083938695 * np.exp(-0.0900900900900901 * v)),
        310000.0 * np.exp(0.3485 * v) + 2.7 * np.exp(0.079 * v),
    )
    th = 1.0 / (ah + bh)
    dh_dt = (-h + hss) / th
    dh_dt_linearized = -1 / th
    states[12] = (-1 + np.exp(dt * dh_dt_linearized)) * dh_dt / dh_dt_linearized + h
    aj = np.where(
        v >= -40,
        0,
        (37.78 + v)
        * (-25428.0 * np.exp(0.2444 * v) - 6.948e-06 * np.exp(-0.04391 * v))
        / (1 + 50262745825.95399 * np.exp(0.311 * v)),
    )
    bj = np.where(
        v >= -40,
        0.6 * np.exp(0.057 * v) / (1 + 0.040762203978366204 * np.exp(-0.1 * v)),
        0.02424
        * np.exp(-0.01052 * v)
        / (1 + 0.003960868339904256 * np.exp(-0.1378 * v)),
    )
    jss = hss
    tj = 1.0 / (aj + bj)
    dj_dt = (-j + jss) / tj
    dj_dt_linearized = -1 / tj
    states[13] = (-1 + np.exp(dt * dj_dt_linearized)) * dj_dt / dj_dt_linearized + j
    hssp = 1.0 / (
        (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
        * (1 + 34112.38799331305 * np.exp(0.13458950201884254 * v))
    )
    dhp_dt = (-hp + hssp) / th
    dhp_dt_linearized = -1 / th
    states[14] = (-1 + np.exp(dt * dhp_dt_linearized)) * dhp_dt / dhp_dt_linearized + hp
    tjp = 1.46 * tj
    djp_dt = (-jp + jss) / tjp
    djp_dt_linearized = -1 / tjp
    states[15] = (-1 + np.exp(dt * djp_dt_linearized)) * djp_dt / djp_dt_linearized + jp
    fINap = 1.0 / (1 + KmCaMK / CaMKa)
    INa = GNa * (m * m * m) * (-ENa + v) * ((1 - fINap) * h * j + fINap * hp * jp)

    # Expressions for the L component
    mLss = 1.0 / (1 + 0.000291579585635531 * np.exp(-0.18996960486322187 * v))
    tmL = 0.1292 * np.exp(
        -(
            (2.9465894465894467 + 0.06435006435006435 * v)
            * (2.9465894465894467 + 0.06435006435006435 * v)
        )
    ) + 0.06487 * np.exp(
        -(
            (-0.09434663536776214 + 0.019561815336463225 * v)
            * (-0.09434663536776214 + 0.019561815336463225 * v)
        )
    )
    dmL_dt = (-mL + mLss) / tmL
    dmL_dt_linearized = -1 / tmL
    states[16] = (-1 + np.exp(dt * dmL_dt_linearized)) * dmL_dt / dmL_dt_linearized + mL
    hLss = 1.0 / (1 + 120578.15595522427 * np.exp(0.13354700854700854 * v))
    dhL_dt = (-hL + hLss) / thL
    dhL_dt_linearized = -1 / thL
    states[17] = (-1 + np.exp(dt * dhL_dt_linearized)) * dhL_dt / dhL_dt_linearized + hL
    hLssp = 1.0 / (1 + 275969.2903869871 * np.exp(0.13354700854700854 * v))
    thLp = 3 * thL
    dhLp_dt = (-hLp + hLssp) / thLp
    dhLp_dt_linearized = -1 / thLp
    states[18] = (
        -1 + np.exp(dt * dhLp_dt_linearized)
    ) * dhLp_dt / dhLp_dt_linearized + hLp
    GNaL = np.where(celltype == 1, 0.6 * GNaL_b, GNaL_b)
    fINaLp = 1.0 / (1 + KmCaMK / CaMKa)
    INaL = (-ENa + v) * ((1 - fINaLp) * hL + fINaLp * hLp) * GNaL * mL

    # Expressions for the Ito component
    ass = 1.0 / (
        1
        + 2.6316508161673635
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    ta = 1.0515 / (
        1.0
        / (
            1.2089
            + 2.2621017070578837
            * np.exp(-0.03403513787634354 * EKshift - 0.03403513787634354 * v)
        )
        + 3.5
        / (
            1
            + 30.069572727397507
            * np.exp(0.03403513787634354 * EKshift + 0.03403513787634354 * v)
        )
    )
    da_dt = (-a + ass) / ta
    da_dt_linearized = -1 / ta
    states[19] = (-1 + np.exp(dt * da_dt_linearized)) * da_dt / da_dt_linearized + a
    iss = 1.0 / (
        1
        + 2194.970764538301
        * np.exp(0.17510068289266328 * EKshift + 0.17510068289266328 * v)
    )
    delta_epi = np.where(
        celltype == 1, 1 - 0.95 / (1 + np.exp(14 + EKshift / 5 + v / 5)), 1
    )
    tiF_b = 4.562 + 1.0 / (
        0.3933 * np.exp(-1 - EKshift / 100 - v / 100)
        + 1.6300896349780942
        * np.exp(0.06027727546714889 * EKshift + 0.06027727546714889 * v)
    )
    tiS_b = 23.62 + 1.0 / (
        0.00027617763953377436
        * np.exp(-0.01693480101608806 * EKshift - 0.01693480101608806 * v)
        + 0.024208962804604526
        * np.exp(0.12377769525931426 * EKshift + 0.12377769525931426 * v)
    )
    tiF = delta_epi * tiF_b
    tiS = delta_epi * tiS_b
    AiF = 1.0 / (
        1
        + 0.24348537187522867
        * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
    )
    AiS = 1 - AiF
    diF_dt = (-iF + iss) / tiF
    diF_dt_linearized = -1 / tiF
    states[20] = (-1 + np.exp(dt * diF_dt_linearized)) * diF_dt / diF_dt_linearized + iF
    diS_dt = (-iS + iss) / tiS
    diS_dt_linearized = -1 / tiS
    states[21] = (-1 + np.exp(dt * diS_dt_linearized)) * diS_dt / diS_dt_linearized + iS
    i = AiF * iF + AiS * iS
    assp = 1.0 / (
        1
        + 5.167428462230666
        * np.exp(-0.06747638326585695 * EKshift - 0.06747638326585695 * v)
    )
    dap_dt = (-ap + assp) / ta
    dap_dt_linearized = -1 / ta
    states[22] = (-1 + np.exp(dt * dap_dt_linearized)) * dap_dt / dap_dt_linearized + ap
    dti_develop = 1.354 + 0.0001 / (
        2.6591269045230603e-05
        * np.exp(0.06293266205160478 * EKshift + 0.06293266205160478 * v)
        + 4.5541779737128264e24
        * np.exp(-4.642525533890436 * EKshift - 4.642525533890436 * v)
    )
    dti_recover = 1 - 0.5 / (1 + 33.11545195869231 * np.exp(0.05 * EKshift + 0.05 * v))
    tiFp = dti_develop * dti_recover * tiF
    tiSp = dti_develop * dti_recover * tiS
    diFp_dt = (-iFp + iss) / tiFp
    diFp_dt_linearized = -1 / tiFp
    states[23] = (
        -1 + np.exp(dt * diFp_dt_linearized)
    ) * diFp_dt / diFp_dt_linearized + iFp
    diSp_dt = (-iSp + iss) / tiSp
    diSp_dt_linearized = -1 / tiSp
    states[24] = (
        -1 + np.exp(dt * diSp_dt_linearized)
    ) * diSp_dt / diSp_dt_linearized + iSp
    ip = AiF * iFp + AiS * iSp
    Gto = np.where(celltype == 1, 2 * Gto_b, np.where(celltype == 2, 2 * Gto_b, Gto_b))
    fItop = 1.0 / (1 + KmCaMK / CaMKa)
    Ito = (-EK + v) * ((1 - fItop) * a * i + ap * fItop * ip) * Gto

    # Expressions for the IKr component
    alpha = 0.1161 * np.exp(0.299 * vfrt)
    beta = 0.2442 * np.exp(-1.604 * vfrt)
    alpha_2 = 0.0578 * np.exp(0.971 * vfrt)
    beta_2 = 0.000349 * np.exp(-1.062 * vfrt)
    alpha_i = 0.2533 * np.exp(0.5953 * vfrt)
    beta_i = 0.06525 * np.exp(-0.8209 * vfrt)
    alpha_C2ToI = 5.2e-05 * np.exp(1.525 * vfrt)
    beta_ItoC2 = alpha_C2ToI * beta_2 * beta_i / (alpha_2 * alpha_i)
    dC3_dt = C2 * beta - C3 * alpha
    dC3_dt_linearized = -alpha
    states[25] = C3 + np.where(
        np.abs(dC3_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dC3_dt_linearized)) * dC3_dt / dC3_dt_linearized,
        dt * dC3_dt,
    )
    dC2_dt = beta_1 * C1 + C3 * alpha - (alpha_1 + beta) * C2
    dC2_dt_linearized = -alpha_1 - beta
    states[26] = C2 + np.where(
        np.abs(dC2_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dC2_dt_linearized)) * dC2_dt / dC2_dt_linearized,
        dt * dC2_dt,
    )
    dC1_dt = (
        alpha_1 * C2
        + I * beta_ItoC2
        + O * beta_2
        - (beta_1 + alpha_2 + alpha_C2ToI) * C1
    )
    dC1_dt_linearized = -beta_1 - alpha_2 - alpha_C2ToI
    states[27] = C1 + np.where(
        np.abs(dC1_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dC1_dt_linearized)) * dC1_dt / dC1_dt_linearized,
        dt * dC1_dt,
    )
    dO_dt = C1 * alpha_2 + I * beta_i - (alpha_i + beta_2) * O
    dO_dt_linearized = -alpha_i - beta_2
    states[28] = O + np.where(
        np.abs(dO_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dO_dt_linearized)) * dO_dt / dO_dt_linearized,
        dt * dO_dt,
    )
    dI_dt = C1 * alpha_C2ToI + O * alpha_i - (beta_ItoC2 + beta_i) * I
    dI_dt_linearized = -beta_ItoC2 - beta_i
    states[29] = I + np.where(
        np.abs(dI_dt_linearized) > 1e-08,
        (-1 + np.exp(dt * dI_dt_linearized)) * dI_dt / dI_dt_linearized,
        dt * dI_dt,
    )
    GKr = np.where(
        celltype == 1, 1.3 * GKr_b, np.where(celltype == 2, 0.8 * GKr_b, GKr_b)
    )
    IKr = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GKr * O / 5

    # Expressions for the IKs component
    xs1ss = 1.0 / (1 + 0.27288596035656526 * np.exp(-0.11195700850873264 * v))
    txs1 = 817.3 + 1.0 / (
        0.003504067763074858 * np.exp(0.056179775280898875 * v)
        + 0.001292 * np.exp(-21 / 23 - v / 230)
    )
    dxs1_dt = (-xs1 + xs1ss) / txs1
    dxs1_dt_linearized = -1 / txs1
    states[30] = (
        -1 + np.exp(dt * dxs1_dt_linearized)
    ) * dxs1_dt / dxs1_dt_linearized + xs1
    xs2ss = xs1ss
    txs2 = 1.0 / (
        0.0022561357010639103 * np.exp(-v / 31) + 0.01 * np.exp(-5 / 2 + v / 20)
    )
    dxs2_dt = (-xs2 + xs2ss) / txs2
    dxs2_dt_linearized = -1 / txs2
    states[31] = (
        -1 + np.exp(dt * dxs2_dt_linearized)
    ) * dxs2_dt / dxs2_dt_linearized + xs2
    KsCa = 1 + 0.6 / (1 + 6.481821026062645e-07 * np.power(1.0 / cai, 1.4))
    GKs = np.where(celltype == 1, 1.4 * GKs_b, GKs_b)
    IKs = (-EKs + v) * GKs * KsCa * xs1 * xs2

    # Expressions for the IK1 component
    aK1 = 4.094 / (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
    bK1 = (
        12.621629724407278 * np.exp(0.0674 * v - 0.0674 * EK)
        + 1.1196358381249121e-16 * np.exp(0.0618 * v - 0.0618 * EK)
    ) / (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
    K1ss = aK1 / (aK1 + bK1)
    GK1 = np.where(
        celltype == 1, 1.2 * GK1_b, np.where(celltype == 2, 1.3 * GK1_b, GK1_b)
    )
    IK1 = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GK1 * K1ss / 5

    # Expressions for the IKb component
    xkb = 1.0 / (1 + 1.57503502085457 * np.exp(-0.04168907454423419 * v))
    GKb = np.where(celltype == 1, 0.6 * GKb_b, GKb_b)
    IKb = (-EK + v) * GKb * xkb

    # Expressions for the ICl component
    IClCa_junc = Fjunc * GClCa * (-EClss + v) / (1 + KdClCa / cass)
    IClCa_sl = GClCa * (1 - Fjunc) * (-ECl + v) / (1 + KdClCa / cai)
    IClCa = IClCa_junc + IClCa_sl
    IClb = GClb * (-ECl + v)

    # Expressions for the ICab component
    ICab = (
        4
        * PCab
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )

    # Expressions for the Ryr component
    a_rel = 0.5 * bt
    Jrel_inf_b = -ICaL_ss * a_rel / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    Jrel_inf = np.where(celltype == 2, 1.7 * Jrel_inf_b, Jrel_inf_b)
    tau_rel_b = bt / (1 + 0.0123 / cajsr)
    tau_rel = np.where(tau_rel_b < 0.001, 0.001, tau_rel_b)
    dJrel_np_dt = (-Jrel_np + Jrel_inf) / tau_rel
    dJrel_np_dt_linearized = -1 / tau_rel
    states[32] = (
        -1 + np.exp(dt * dJrel_np_dt_linearized)
    ) * dJrel_np_dt / dJrel_np_dt_linearized + Jrel_np
    btp = 1.25 * bt
    a_relp = 0.5 * btp
    Jrel_infp_b = -ICaL_ss * a_relp / (1 + np.power(cajsr_half, 8) / np.power(cajsr, 8))
    Jrel_infp = np.where(celltype == 2, 1.7 * Jrel_infp_b, Jrel_infp_b)
    tau_relp_b = btp / (1 + 0.0123 / cajsr)
    tau_relp = np.where(tau_relp_b < 0.001, 0.001, tau_relp_b)
    dJrel_p_dt = (-Jrel_p + Jrel_infp) / tau_relp
    dJrel_p_dt_linearized = -1 / tau_relp
    states[33] = (
        -1 + np.exp(dt * dJrel_p_dt_linearized)
    ) * dJrel_p_dt / dJrel_p_dt_linearized + Jrel_p
    fJrelp = 1.0 / (1 + KmCaMK / CaMKa)
    Jrel = Jrel_b * ((1 - fJrelp) * Jrel_np + Jrel_p * fJrelp)

    # Expressions for the Membrane component
    Istim = np.where(
        np.logical_and(
            t >= i_Stim_Start,
            t
            - i_Stim_Start
            - i_Stim_Period * np.floor((t - i_Stim_Start) / i_Stim_Period)
            <= i_Stim_PulseDuration,
        ),
        i_Stim_Amplitude,
        0,
    )
    dv_dt = (
        -ICaK
        - ICaL
        - ICaNa
        - ICab
        - IClCa
        - IClb
        - IK1
        - IKb
        - IKr
        - IKs
        - INa
        - INaCa_i
        - INaCa_ss
        - INaK
        - INaL
        - INab
        - I_katp
        - IpCa
        - Istim
        - Ito
    )
    dAfcaf_dv = (
        -0.06
        * np.exp(-1 + v / 10)
        / ((1 + np.exp(-1 + v / 10)) * (1 + np.exp(-1 + v / 10)))
    )
    dAiF_dv = (
        -0.0016103529885927823
        * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
        / (
            (
                1
                + 0.24348537187522867
                * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
            )
            * (
                1
                + 0.24348537187522867
                * np.exp(0.006613756613756614 * EKshift + 0.006613756613756614 * v)
            )
        )
    )
    dE1_dx1 = 1.0 / (x1 + x2 + x3 + x4) - x1 / (
        (x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4)
    )
    dE1_dx2 = -x1 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE1_dx3 = -x1 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE1_dx4 = -x1 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE1_i_dx1_i = 1.0 / (x1_i + x2_i + x3_i + x4_i) - x1_i / (
        (x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i)
    )
    dE1_i_dx2_i = -x1_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE1_i_dx3_i = -x1_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE1_i_dx4_i = -x1_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE1_ss_dx1_ss = 1.0 / (x1_ss + x2_ss + x3_ss + x4_ss) - x1_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE1_ss_dx2_ss = -x1_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE1_ss_dx3_ss = -x1_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE1_ss_dx4_ss = -x1_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE2_dx1 = -x2 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE2_dx2 = 1.0 / (x1 + x2 + x3 + x4) - x2 / (
        (x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4)
    )
    dE2_dx3 = -x2 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE2_dx4 = -x2 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE2_i_dx1_i = -x2_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE2_i_dx2_i = 1.0 / (x1_i + x2_i + x3_i + x4_i) - x2_i / (
        (x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i)
    )
    dE2_i_dx3_i = -x2_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE2_i_dx4_i = -x2_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE2_ss_dx1_ss = -x2_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE2_ss_dx2_ss = 1.0 / (x1_ss + x2_ss + x3_ss + x4_ss) - x2_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE2_ss_dx3_ss = -x2_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE2_ss_dx4_ss = -x2_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE3_dx1 = -x3 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE3_dx2 = -x3 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE3_dx3 = 1.0 / (x1 + x2 + x3 + x4) - x3 / (
        (x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4)
    )
    dE3_dx4 = -x3 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE3_i_dx1_i = -x3_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE3_i_dx2_i = -x3_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE3_i_dx3_i = 1.0 / (x1_i + x2_i + x3_i + x4_i) - x3_i / (
        (x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i)
    )
    dE3_i_dx4_i = -x3_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE3_ss_dx1_ss = -x3_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE3_ss_dx2_ss = -x3_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE3_ss_dx3_ss = 1.0 / (x1_ss + x2_ss + x3_ss + x4_ss) - x3_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE3_ss_dx4_ss = -x3_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE4_dx1 = -x4 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE4_dx2 = -x4 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE4_dx3 = -x4 / ((x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4))
    dE4_dx4 = 1.0 / (x1 + x2 + x3 + x4) - x4 / (
        (x1 + x2 + x3 + x4) * (x1 + x2 + x3 + x4)
    )
    dE4_i_dx1_i = -x4_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE4_i_dx2_i = -x4_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE4_i_dx3_i = -x4_i / ((x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i))
    dE4_i_dx4_i = 1.0 / (x1_i + x2_i + x3_i + x4_i) - x4_i / (
        (x1_i + x2_i + x3_i + x4_i) * (x1_i + x2_i + x3_i + x4_i)
    )
    dE4_ss_dx1_ss = -x4_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE4_ss_dx2_ss = -x4_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE4_ss_dx3_ss = -x4_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dE4_ss_dx4_ss = 1.0 / (x1_ss + x2_ss + x3_ss + x4_ss) - x4_ss / (
        (x1_ss + x2_ss + x3_ss + x4_ss) * (x1_ss + x2_ss + x3_ss + x4_ss)
    )
    dICaK_i_dPhiCaK_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaK * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaKp * d * fICaLp
    )
    dICaK_i_dfca = (
        (1 - ICaL_fractionSS) * (1 - fICaLp) * PCaK * PhiCaK_i * d * jca * nca_i
    )
    dICaK_i_dfcap = (1 - ICaL_fractionSS) * PCaKp * PhiCaK_i * d * fICaLp * jca * nca_i
    dICaK_ss_dPhiCaK_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaK * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaKp * d * fICaLp
    )
    dICaK_ss_dfca = ICaL_fractionSS * (1 - fICaLp) * PCaK * PhiCaK_ss * d * jca * nca_ss
    dICaK_ss_dfcap = ICaL_fractionSS * PCaKp * PhiCaK_ss * d * fICaLp * jca * nca_ss
    dICaL_i_dPhiCaL_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCa * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCap * d * fICaLp
    )
    dICaL_i_dfca = (
        (1 - ICaL_fractionSS) * (1 - fICaLp) * PCa * PhiCaL_i * d * jca * nca_i
    )
    dICaL_i_dfcap = (1 - ICaL_fractionSS) * PCap * PhiCaL_i * d * fICaLp * jca * nca_i
    dICaL_ss_dPhiCaL_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCa * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCap * d * fICaLp
    )
    dICaL_ss_dfca = ICaL_fractionSS * (1 - fICaLp) * PCa * PhiCaL_ss * d * jca * nca_ss
    dICaL_ss_dfcap = ICaL_fractionSS * PCap * PhiCaL_ss * d * fICaLp * jca * nca_ss
    dICaNa_i_dPhiCaNa_i = (1 - ICaL_fractionSS) * (
        (1 - fICaLp) * ((1 - nca_i) * f + fca * jca * nca_i) * PCaNa * d
        + ((1 - nca_i) * fp + fcap * jca * nca_i) * PCaNap * d * fICaLp
    )
    dICaNa_i_dfca = (
        (1 - ICaL_fractionSS) * (1 - fICaLp) * PCaNa * PhiCaNa_i * d * jca * nca_i
    )
    dICaNa_i_dfcap = (
        (1 - ICaL_fractionSS) * PCaNap * PhiCaNa_i * d * fICaLp * jca * nca_i
    )
    dICaNa_ss_dPhiCaNa_ss = ICaL_fractionSS * (
        (1 - fICaLp) * ((1 - nca_ss) * f + fca * jca * nca_ss) * PCaNa * d
        + ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCaNap * d * fICaLp
    )
    dICaNa_ss_dfca = (
        ICaL_fractionSS * (1 - fICaLp) * PCaNa * PhiCaNa_ss * d * jca * nca_ss
    )
    dICaNa_ss_dfcap = ICaL_fractionSS * PCaNap * PhiCaNa_ss * d * fICaLp * jca * nca_ss
    dICab_dvffrt = (
        4
        * PCab
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        / (-1 + np.exp(2 * vfrt))
    )
    dICab_dvfrt = -8 * PCab * (
        -cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai
    ) * np.exp(2 * vfrt) * vffrt / (
        (-1 + np.exp(2 * vfrt)) * (-1 + np.exp(2 * vfrt))
    ) + 8 * PCab * cai * np.exp(
        2 * vfrt
    ) * gamma_cai * vffrt / (
        -1 + np.exp(2 * vfrt)
    )
    dIClCa_junc_dv = Fjunc * GClCa / (1 + KdClCa / cass)
    dIClCa_sl_dv = GClCa * (1 - Fjunc) / (1 + KdClCa / cai)
    dIK1_dK1ss = np.sqrt(5) * np.sqrt(ko) * (-EK + v) * GK1 / 5
    dK1ss_daK1 = 1.0 / (aK1 + bK1) - aK1 / ((aK1 + bK1) * (aK1 + bK1))
    dK1ss_dbK1 = -aK1 / ((aK1 + bK1) * (aK1 + bK1))
    daK1_dv = (
        -0.0011435228161993972
        * np.exp(0.1217 * v - 0.1217 * EK)
        / (
            (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
            * (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
        )
    )
    dbK1_dv = (
        6.919349479611957e-18 * np.exp(0.0618 * v - 0.0618 * EK)
        + 0.8506978434250505 * np.exp(0.0674 * v - 0.0674 * EK)
    ) / (
        1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v)
    ) + 0.016099950792318272 * (
        12.621629724407278 * np.exp(0.0674 * v - 0.0674 * EK)
        + 1.1196358381249121e-16 * np.exp(0.0618 * v - 0.0618 * EK)
    ) * np.exp(
        0.1629 * EK - 0.1629 * v
    ) / (
        (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
        * (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
    )
    dIK1_dv = (
        np.sqrt(5) * np.sqrt(ko) * GK1 * K1ss / 5
        + np.sqrt(5)
        * np.sqrt(ko)
        * (-EK + v)
        * (dK1ss_daK1 * daK1_dv + dK1ss_dbK1 * dbK1_dv)
        * GK1
        / 5
    )
    dxkb_dv = (
        0.06566175239418562
        * np.exp(-0.04168907454423419 * v)
        / (
            (1 + 1.57503502085457 * np.exp(-0.04168907454423419 * v))
            * (1 + 1.57503502085457 * np.exp(-0.04168907454423419 * v))
        )
    )
    dIKb_dv = GKb * xkb + (-EK + v) * GKb * dxkb_dv
    dIKb_dxkb = (-EK + v) * GKb
    dIKr_dv = np.sqrt(5) * np.sqrt(ko) * GKr * O / 5
    dIKs_dv = GKs * KsCa * xs1 * xs2
    dINa_dv = GNa * (m * m * m) * ((1 - fINap) * h * j + fINap * hp * jp)
    dINaCa_i_dJncxCa_i = zca * (1 - INaCa_fractionSS) * Gncx * allo_i
    dINaCa_i_dJncxNa_i = zna * (1 - INaCa_fractionSS) * Gncx * allo_i
    dINaCa_ss_dJncxCa_ss = INaCa_fractionSS * zca * Gncx * allo_ss
    dINaCa_ss_dJncxNa_ss = INaCa_fractionSS * zna * Gncx * allo_ss
    dINaL_dv = ((1 - fINaLp) * hL + fINaLp * hLp) * GNaL * mL
    dINab_dvffrt = PNab * (-nao + np.exp(vfrt) * nai) / (-1 + np.exp(vfrt))
    dINab_dvfrt = PNab * np.exp(vfrt) * nai * vffrt / (-1 + np.exp(vfrt)) - PNab * (
        -nao + np.exp(vfrt) * nai
    ) * np.exp(vfrt) * vffrt / ((-1 + np.exp(vfrt)) * (-1 + np.exp(vfrt)))
    dI_katp_dv = fkatp * gkatp * akik * bkik
    dIto_di = (1 - fItop) * (-EK + v) * Gto * a
    dIto_dip = (-EK + v) * Gto * ap * fItop
    di_dAiF = -iS + iF
    dip_dAiF = -iSp + iFp
    dIto_dv = ((1 - fItop) * a * i + ap * fItop * ip) * Gto + (-EK + v) * (
        (1 - fItop) * (dAiF_dv * di_dAiF - dAiF_dv * iS) * a
        + (dAiF_dv * dip_dAiF - dAiF_dv * iSp) * ap * fItop
    ) * Gto
    dx1_da1 = a2 * a4 + a2 * b3
    dx4_da1 = a3 * a4 + a4 * b2 + b2 * b3
    dJnakK_da1 = (
        -2 * E3
        - 2 * (dE3_dx1 * dx1_da1 + dE3_dx4 * dx4_da1 + a2 * a3 * dE3_dx2) * a1
        + 2 * (dE4_dx1 * dx1_da1 + dE4_dx4 * dx4_da1 + a2 * a3 * dE4_dx2) * b1
    )
    dx2_da3 = a1 * a2 + a2 * b4 + b1 * b4
    dx3_da3 = a2 * a4 + a4 * b1
    dJnakNa_da3 = (
        3 * E1
        - 3 * (dE2_dx2 * dx2_da3 + dE2_dx3 * dx3_da3 + a1 * a4 * dE2_dx4) * b3
        + 3 * (dE1_dx2 * dx2_da3 + dE1_dx3 * dx3_da3 + a1 * a4 * dE1_dx4) * a3
    )
    dx3_i_dk3_i = (k6_i + k7_i) * k1_i + k6_i * k8_i
    dx4_i_dk3_i = (k1_i + k8_i) * k5_i
    dJncxNa_i_dk3pp_i = (
        -E2_i
        + (
            dE3_i_dx3_i * dx3_i_dk3_i
            + dE3_i_dx4_i * dx4_i_dk3_i
            + dE3_i_dx1_i * k5_i * k7_i
        )
        * k4pp_i
        - (
            dE2_i_dx3_i * dx3_i_dk3_i
            + dE2_i_dx4_i * dx4_i_dk3_i
            + dE2_i_dx1_i * k5_i * k7_i
        )
        * k3pp_i
        - 3
        * (
            dE1_i_dx3_i * dx3_i_dk3_i
            + dE1_i_dx4_i * dx4_i_dk3_i
            + dE1_i_dx1_i * k5_i * k7_i
        )
        * k8_i
        + 3
        * (
            dE4_i_dx3_i * dx3_i_dk3_i
            + dE4_i_dx4_i * dx4_i_dk3_i
            + dE4_i_dx1_i * k5_i * k7_i
        )
        * k7_i
    )
    dx1_i_dk4_i = (k6_i + k7_i) * k2_i
    dx2_i_dk4_i = (k1_i + k8_i) * k6_i + k1_i * k7_i
    dJncxNa_i_dk4pp_i = (
        (
            dE3_i_dx1_i * dx1_i_dk4_i
            + dE3_i_dx2_i * dx2_i_dk4_i
            + dE3_i_dx4_i * k2_i * k8_i
        )
        * k4pp_i
        - (
            dE2_i_dx1_i * dx1_i_dk4_i
            + dE2_i_dx2_i * dx2_i_dk4_i
            + dE2_i_dx4_i * k2_i * k8_i
        )
        * k3pp_i
        - 3
        * (
            dE1_i_dx1_i * dx1_i_dk4_i
            + dE1_i_dx2_i * dx2_i_dk4_i
            + dE1_i_dx4_i * k2_i * k8_i
        )
        * k8_i
        + 3
        * (
            dE4_i_dx1_i * dx1_i_dk4_i
            + dE4_i_dx2_i * dx2_i_dk4_i
            + dE4_i_dx4_i * k2_i * k8_i
        )
        * k7_i
        + E3_i
    )
    dx1_i_dk7_i = (k2_i + k3_i) * k5_i + k2_i * k4_i
    dx2_i_dk7_i = (k4_i + k5_i) * k1_i
    dJncxNa_i_dk7_i = (
        3 * E4_i
        + (
            dE3_i_dx1_i * dx1_i_dk7_i
            + dE3_i_dx2_i * dx2_i_dk7_i
            + dE3_i_dx3_i * k1_i * k3_i
        )
        * k4pp_i
        - (
            dE2_i_dx1_i * dx1_i_dk7_i
            + dE2_i_dx2_i * dx2_i_dk7_i
            + dE2_i_dx3_i * k1_i * k3_i
        )
        * k3pp_i
        - 3
        * (
            dE1_i_dx1_i * dx1_i_dk7_i
            + dE1_i_dx2_i * dx2_i_dk7_i
            + dE1_i_dx3_i * k1_i * k3_i
        )
        * k8_i
        + 3
        * (
            dE4_i_dx1_i * dx1_i_dk7_i
            + dE4_i_dx2_i * dx2_i_dk7_i
            + dE4_i_dx3_i * k1_i * k3_i
        )
        * k7_i
    )
    dx3_i_dk8_i = (k2_i + k3_i) * k6_i
    dx4_i_dk8_i = (k4_i + k5_i) * k2_i + k3_i * k5_i
    dJncxNa_i_dk8_i = (
        -3 * E1_i
        + (
            dE3_i_dx3_i * dx3_i_dk8_i
            + dE3_i_dx4_i * dx4_i_dk8_i
            + dE3_i_dx2_i * k4_i * k6_i
        )
        * k4pp_i
        - (
            dE2_i_dx3_i * dx3_i_dk8_i
            + dE2_i_dx4_i * dx4_i_dk8_i
            + dE2_i_dx2_i * k4_i * k6_i
        )
        * k3pp_i
        - 3
        * (
            dE1_i_dx3_i * dx3_i_dk8_i
            + dE1_i_dx4_i * dx4_i_dk8_i
            + dE1_i_dx2_i * k4_i * k6_i
        )
        * k8_i
        + 3
        * (
            dE4_i_dx3_i * dx3_i_dk8_i
            + dE4_i_dx4_i * dx4_i_dk8_i
            + dE4_i_dx2_i * k4_i * k6_i
        )
        * k7_i
    )
    dx3_ss_dk3_ss = (k6_ss + k7_ss) * k1_ss + k6_ss * k8_ss
    dx4_ss_dk3_ss = (k1_ss + k8_ss) * k5_ss
    dJncxNa_ss_dk3pp_ss = (
        -E2_ss
        + (
            dE3_ss_dx3_ss * dx3_ss_dk3_ss
            + dE3_ss_dx4_ss * dx4_ss_dk3_ss
            + dE3_ss_dx1_ss * k5_ss * k7_ss
        )
        * k4pp_ss
        - (
            dE2_ss_dx3_ss * dx3_ss_dk3_ss
            + dE2_ss_dx4_ss * dx4_ss_dk3_ss
            + dE2_ss_dx1_ss * k5_ss * k7_ss
        )
        * k3pp_ss
        - 3
        * (
            dE1_ss_dx3_ss * dx3_ss_dk3_ss
            + dE1_ss_dx4_ss * dx4_ss_dk3_ss
            + dE1_ss_dx1_ss * k5_ss * k7_ss
        )
        * k8_ss
        + 3
        * (
            dE4_ss_dx3_ss * dx3_ss_dk3_ss
            + dE4_ss_dx4_ss * dx4_ss_dk3_ss
            + dE4_ss_dx1_ss * k5_ss * k7_ss
        )
        * k7_ss
    )
    dx1_ss_dk4_ss = (k6_ss + k7_ss) * k2_ss
    dx2_ss_dk4_ss = (k1_ss + k8_ss) * k6_ss + k1_ss * k7_ss
    dJncxNa_ss_dk4pp_ss = (
        (
            dE3_ss_dx1_ss * dx1_ss_dk4_ss
            + dE3_ss_dx2_ss * dx2_ss_dk4_ss
            + dE3_ss_dx4_ss * k2_ss * k8_ss
        )
        * k4pp_ss
        - (
            dE2_ss_dx1_ss * dx1_ss_dk4_ss
            + dE2_ss_dx2_ss * dx2_ss_dk4_ss
            + dE2_ss_dx4_ss * k2_ss * k8_ss
        )
        * k3pp_ss
        - 3
        * (
            dE1_ss_dx1_ss * dx1_ss_dk4_ss
            + dE1_ss_dx2_ss * dx2_ss_dk4_ss
            + dE1_ss_dx4_ss * k2_ss * k8_ss
        )
        * k8_ss
        + 3
        * (
            dE4_ss_dx1_ss * dx1_ss_dk4_ss
            + dE4_ss_dx2_ss * dx2_ss_dk4_ss
            + dE4_ss_dx4_ss * k2_ss * k8_ss
        )
        * k7_ss
        + E3_ss
    )
    dx1_ss_dk7_ss = (k2_ss + k3_ss) * k5_ss + k2_ss * k4_ss
    dx2_ss_dk7_ss = (k4_ss + k5_ss) * k1_ss
    dJncxNa_ss_dk7_ss = (
        3 * E4_ss
        + (
            dE3_ss_dx1_ss * dx1_ss_dk7_ss
            + dE3_ss_dx2_ss * dx2_ss_dk7_ss
            + dE3_ss_dx3_ss * k1_ss * k3_ss
        )
        * k4pp_ss
        - (
            dE2_ss_dx1_ss * dx1_ss_dk7_ss
            + dE2_ss_dx2_ss * dx2_ss_dk7_ss
            + dE2_ss_dx3_ss * k1_ss * k3_ss
        )
        * k3pp_ss
        - 3
        * (
            dE1_ss_dx1_ss * dx1_ss_dk7_ss
            + dE1_ss_dx2_ss * dx2_ss_dk7_ss
            + dE1_ss_dx3_ss * k1_ss * k3_ss
        )
        * k8_ss
        + 3
        * (
            dE4_ss_dx1_ss * dx1_ss_dk7_ss
            + dE4_ss_dx2_ss * dx2_ss_dk7_ss
            + dE4_ss_dx3_ss * k1_ss * k3_ss
        )
        * k7_ss
    )
    dx3_ss_dk8_ss = (k2_ss + k3_ss) * k6_ss
    dx4_ss_dk8_ss = (k4_ss + k5_ss) * k2_ss + k3_ss * k5_ss
    dJncxNa_ss_dk8_ss = (
        -3 * E1_ss
        + (
            dE3_ss_dx3_ss * dx3_ss_dk8_ss
            + dE3_ss_dx4_ss * dx4_ss_dk8_ss
            + dE3_ss_dx2_ss * k4_ss * k6_ss
        )
        * k4pp_ss
        - (
            dE2_ss_dx3_ss * dx3_ss_dk8_ss
            + dE2_ss_dx4_ss * dx4_ss_dk8_ss
            + dE2_ss_dx2_ss * k4_ss * k6_ss
        )
        * k3pp_ss
        - 3
        * (
            dE1_ss_dx3_ss * dx3_ss_dk8_ss
            + dE1_ss_dx4_ss * dx4_ss_dk8_ss
            + dE1_ss_dx2_ss * k4_ss * k6_ss
        )
        * k8_ss
        + 3
        * (
            dE4_ss_dx3_ss * dx3_ss_dk8_ss
            + dE4_ss_dx4_ss * dx4_ss_dk8_ss
            + dE4_ss_dx2_ss * k4_ss * k6_ss
        )
        * k7_ss
    )
    dKnai_dvfrt = Knai0 * delta * np.exp(delta * vfrt / 3) / 3
    dKnao_dvfrt = Knao0 * (1 / 3 - delta / 3) * np.exp((1 - delta) * vfrt / 3)
    dPhiCaK_i_dvffrt = (-ko * gamma_ko + np.exp(vfrt) * gamma_ki * ki) / (
        -1 + np.exp(vfrt)
    )
    dPhiCaK_i_dvfrt = -(-ko * gamma_ko + np.exp(vfrt) * gamma_ki * ki) * np.exp(
        vfrt
    ) * vffrt / ((-1 + np.exp(vfrt)) * (-1 + np.exp(vfrt))) + np.exp(
        vfrt
    ) * gamma_ki * ki * vffrt / (
        -1 + np.exp(vfrt)
    )
    dPhiCaK_ss_dvffrt = (-ko * gamma_ko + np.exp(vfrt) * gamma_kss * kss) / (
        -1 + np.exp(vfrt)
    )
    dPhiCaK_ss_dvfrt = -(-ko * gamma_ko + np.exp(vfrt) * gamma_kss * kss) * np.exp(
        vfrt
    ) * vffrt / ((-1 + np.exp(vfrt)) * (-1 + np.exp(vfrt))) + np.exp(
        vfrt
    ) * gamma_kss * kss * vffrt / (
        -1 + np.exp(vfrt)
    )
    dPhiCaL_i_dvffrt = (
        4
        * (-cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai)
        / (-1 + np.exp(2 * vfrt))
    )
    dPhiCaL_i_dvfrt = -8 * (
        -cao * gamma_cao + cai * np.exp(2 * vfrt) * gamma_cai
    ) * np.exp(2 * vfrt) * vffrt / (
        (-1 + np.exp(2 * vfrt)) * (-1 + np.exp(2 * vfrt))
    ) + 8 * cai * np.exp(
        2 * vfrt
    ) * gamma_cai * vffrt / (
        -1 + np.exp(2 * vfrt)
    )
    dPhiCaL_ss_dvffrt = (
        4
        * (-cao * gamma_cao + cass * np.exp(2 * vfrt) * gamma_cass)
        / (-1 + np.exp(2 * vfrt))
    )
    dPhiCaL_ss_dvfrt = -8 * (
        -cao * gamma_cao + cass * np.exp(2 * vfrt) * gamma_cass
    ) * np.exp(2 * vfrt) * vffrt / (
        (-1 + np.exp(2 * vfrt)) * (-1 + np.exp(2 * vfrt))
    ) + 8 * cass * np.exp(
        2 * vfrt
    ) * gamma_cass * vffrt / (
        -1 + np.exp(2 * vfrt)
    )
    dPhiCaNa_i_dvffrt = (-nao * gamma_nao + np.exp(vfrt) * gamma_nai * nai) / (
        -1 + np.exp(vfrt)
    )
    dPhiCaNa_i_dvfrt = -(-nao * gamma_nao + np.exp(vfrt) * gamma_nai * nai) * np.exp(
        vfrt
    ) * vffrt / ((-1 + np.exp(vfrt)) * (-1 + np.exp(vfrt))) + np.exp(
        vfrt
    ) * gamma_nai * nai * vffrt / (
        -1 + np.exp(vfrt)
    )
    dPhiCaNa_ss_dvffrt = (-nao * gamma_nao + np.exp(vfrt) * gamma_nass * nass) / (
        -1 + np.exp(vfrt)
    )
    dPhiCaNa_ss_dvfrt = -(-nao * gamma_nao + np.exp(vfrt) * gamma_nass * nass) * np.exp(
        vfrt
    ) * vffrt / ((-1 + np.exp(vfrt)) * (-1 + np.exp(vfrt))) + np.exp(
        vfrt
    ) * gamma_nass * nass * vffrt / (
        -1 + np.exp(vfrt)
    )
    da1_dKnai = -3 * k1p * (nai * nai * nai) / (
        (
            -1
            + ((1 + ki / Kki) * (1 + ki / Kki))
            + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
        )
        * np.power(Knai, 4)
    ) + 3 * k1p * ((1 + nai / Knai) * (1 + nai / Knai)) * np.power(nai, 4) / (
        (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
        )
        * np.power(Knai, 5)
    )
    da3_dKnao = (
        3
        * k3p
        * nao
        * (ko * ko)
        * ((1 + nao / Knao) * (1 + nao / Knao))
        / (
            (Kko * Kko)
            * (
                (
                    -1
                    + ((1 + ko / Kko) * (1 + ko / Kko))
                    + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
                )
                * (
                    -1
                    + ((1 + ko / Kko) * (1 + ko / Kko))
                    + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
                )
            )
            * (Knao * Knao)
        )
    )
    db2_dKnao = -3 * k2m * (nao * nao * nao) / (
        (
            -1
            + ((1 + ko / Kko) * (1 + ko / Kko))
            + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
        )
        * np.power(Knao, 4)
    ) + 3 * k2m * np.power(nao, 4) * ((1 + nao / Knao) * (1 + nao / Knao)) / (
        (
            (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
            * (
                -1
                + ((1 + ko / Kko) * (1 + ko / Kko))
                + ((1 + nao / Knao) * (1 + nao / Knao) * (1 + nao / Knao))
            )
        )
        * np.power(Knao, 5)
    )
    db4_dKnai = (
        3
        * k4m
        * ((1 + nai / Knai) * (1 + nai / Knai))
        * (ki * ki)
        * nai
        / (
            (Kki * Kki)
            * (
                (
                    -1
                    + ((1 + ki / Kki) * (1 + ki / Kki))
                    + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
                )
                * (
                    -1
                    + ((1 + ki / Kki) * (1 + ki / Kki))
                    + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
                )
            )
            * (Knai * Knai)
        )
    )
    dfca_dAfcaf = -fcas + fcaf
    dfcap_dAfcaf = -fcas + fcafp
    dh1_i_dhna = nai / kna3
    dh1_ss_dhna = nass / kna3
    dh2_i_dh1_i = -hna * nai / (kna3 * (h1_i * h1_i))
    dh2_i_dhna = nai / (kna3 * h1_i) - dh1_i_dhna * hna * nai / (kna3 * (h1_i * h1_i))
    dh2_ss_dh1_ss = -hna * nass / (kna3 * (h1_ss * h1_ss))
    dh2_ss_dhna = nass / (kna3 * h1_ss) - dh1_ss_dhna * hna * nass / (
        kna3 * (h1_ss * h1_ss)
    )
    dh3_i_dh1_i = -1 / (h1_i * h1_i)
    dh3_ss_dh1_ss = -1 / (h1_ss * h1_ss)
    dh7_i_dhna = -nao / (kna3 * (hna * hna))
    dh7_ss_dhna = -nao / (kna3 * (hna * hna))
    dh8_i_dh7_i = -nao / (kna3 * (h7_i * h7_i) * hna)
    dh8_i_dhna = -nao / (kna3 * h7_i * (hna * hna)) - nao * dh7_i_dhna / (
        kna3 * (h7_i * h7_i) * hna
    )
    dh8_ss_dh7_ss = -nao / (kna3 * (h7_ss * h7_ss) * hna)
    dh8_ss_dhna = -nao / (kna3 * h7_ss * (hna * hna)) - nao * dh7_ss_dhna / (
        kna3 * (h7_ss * h7_ss) * hna
    )
    dh9_i_dh7_i = -1 / (h7_i * h7_i)
    dh9_ss_dh7_ss = -1 / (h7_ss * h7_ss)
    dhca_dvfrt = qca * np.exp(qca * vfrt)
    dhna_dvfrt = qna * np.exp(qna * vfrt)
    dk4p_i_dh3_i = wca / hca
    dk4p_i_dhca = -wca * h3_i / (hca * hca)
    dk4p_ss_dh3_ss = wca / hca
    dk4p_ss_dhca = -wca * h3_ss / (hca * hca)
    dvffrt_dv = (F * F) / (R * T)
    dvfrt_dv = F / (R * T)
    dx1_db4 = a2 * b3 + b2 * b3
    dx2_db4 = a2 * a3 + a3 * b1 + b1 * b2
    dx3_db2 = a4 * b1 + b1 * b3
    dx4_db2 = a1 * a4 + a1 * b3 + b3 * b4
    dv_dt_linearized = (
        -GClb
        - dIClCa_junc_dv
        - dIClCa_sl_dv
        - dIK1_dv
        - dIKb_dv
        - dIKr_dv
        - dIKs_dv
        - dINaL_dv
        - dINa_dv
        - dI_katp_dv
        - dIto_dv
        - (
            (
                (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk3_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk8_i
                    * h11_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_i
                    * k1_i
                    * k3_i
                )
                * dE2_i_dx3_i
                + (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk3_i
                    + (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * k2_i
                    * k8_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk8_i
                    * h11_i
                )
                * dE2_i_dx4_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx1_i_dk4_i
                    + (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_i
                    * k7_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_i_dk7_i
                    * h5_i
                )
                * dE2_i_dx1_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx2_i_dk4_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_i_dk7_i
                    * h5_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_i
                    * k4_i
                    * k6_i
                )
                * dE2_i_dx2_i
            )
            * k2_i
            - (
                (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk3_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk8_i
                    * h11_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_i
                    * k1_i
                    * k3_i
                )
                * dE1_i_dx3_i
                + (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk3_i
                    + (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * k2_i
                    * k8_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk8_i
                    * h11_i
                )
                * dE1_i_dx4_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx1_i_dk4_i
                    + (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_i
                    * k7_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_i_dk7_i
                    * h5_i
                )
                * dE1_i_dx1_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx2_i_dk4_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_i_dk7_i
                    * h5_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_i
                    * k4_i
                    * k6_i
                )
                * dE1_i_dx2_i
            )
            * k1_i
        )
        * dINaCa_i_dJncxCa_i
        - (
            (
                (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk3_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk8_ss
                    * h11_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_ss
                    * k1_ss
                    * k3_ss
                )
                * dE2_ss_dx3_ss
                + (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk3_ss
                    + (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * k2_ss
                    * k8_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk8_ss
                    * h11_ss
                )
                * dE2_ss_dx4_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_ss
                    * k7_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_ss_dk7_ss
                    * h5_ss
                )
                * dE2_ss_dx1_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx2_ss_dk4_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_ss_dk7_ss
                    * h5_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_ss
                    * k4_ss
                    * k6_ss
                )
                * dE2_ss_dx2_ss
            )
            * k2_ss
            - (
                (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk3_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk8_ss
                    * h11_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_ss
                    * k1_ss
                    * k3_ss
                )
                * dE1_ss_dx3_ss
                + (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk3_ss
                    + (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * k2_ss
                    * k8_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk8_ss
                    * h11_ss
                )
                * dE1_ss_dx4_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_ss
                    * k7_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_ss_dk7_ss
                    * h5_ss
                )
                * dE1_ss_dx1_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx2_ss_dk4_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_ss_dk7_ss
                    * h5_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_ss
                    * k4_ss
                    * k6_ss
                )
                * dE1_ss_dx2_ss
            )
            * k1_ss
        )
        * dINaCa_ss_dJncxCa_ss
        - (dAfcaf_dv * dfca_dAfcaf - dAfcaf_dv * fcas) * dICaK_i_dfca
        - (dAfcaf_dv * dfca_dAfcaf - dAfcaf_dv * fcas) * dICaK_ss_dfca
        - (dAfcaf_dv * dfca_dAfcaf - dAfcaf_dv * fcas) * dICaL_i_dfca
        - (dAfcaf_dv * dfca_dAfcaf - dAfcaf_dv * fcas) * dICaL_ss_dfca
        - (dAfcaf_dv * dfca_dAfcaf - dAfcaf_dv * fcas) * dICaNa_i_dfca
        - (dAfcaf_dv * dfca_dAfcaf - dAfcaf_dv * fcas) * dICaNa_ss_dfca
        - (dAfcaf_dv * dfcap_dAfcaf - dAfcaf_dv * fcas) * dICaK_i_dfcap
        - (dAfcaf_dv * dfcap_dAfcaf - dAfcaf_dv * fcas) * dICaK_ss_dfcap
        - (dAfcaf_dv * dfcap_dAfcaf - dAfcaf_dv * fcas) * dICaL_i_dfcap
        - (dAfcaf_dv * dfcap_dAfcaf - dAfcaf_dv * fcas) * dICaL_ss_dfcap
        - (dAfcaf_dv * dfcap_dAfcaf - dAfcaf_dv * fcas) * dICaNa_i_dfcap
        - (dAfcaf_dv * dfcap_dAfcaf - dAfcaf_dv * fcas) * dICaNa_ss_dfcap
        - (dAiF_dv * di_dAiF - dAiF_dv * iS) * dIto_di
        - (dAiF_dv * dip_dAiF - dAiF_dv * iSp) * dIto_dip
        - (dK1ss_daK1 * daK1_dv + dK1ss_dbK1 * dbK1_dv) * dIK1_dK1ss
        - (dPhiCaK_i_dvffrt * dvffrt_dv + dPhiCaK_i_dvfrt * dvfrt_dv)
        * dICaK_i_dPhiCaK_i
        - (dPhiCaK_ss_dvffrt * dvffrt_dv + dPhiCaK_ss_dvfrt * dvfrt_dv)
        * dICaK_ss_dPhiCaK_ss
        - (dPhiCaL_i_dvffrt * dvffrt_dv + dPhiCaL_i_dvfrt * dvfrt_dv)
        * dICaL_i_dPhiCaL_i
        - (dPhiCaL_ss_dvffrt * dvffrt_dv + dPhiCaL_ss_dvfrt * dvfrt_dv)
        * dICaL_ss_dPhiCaL_ss
        - (dPhiCaNa_i_dvffrt * dvffrt_dv + dPhiCaNa_i_dvfrt * dvfrt_dv)
        * dICaNa_i_dPhiCaNa_i
        - (dPhiCaNa_ss_dvffrt * dvffrt_dv + dPhiCaNa_ss_dvfrt * dvfrt_dv)
        * dICaNa_ss_dPhiCaNa_ss
        - (
            (
                (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk3_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk8_i
                    * h11_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_i
                    * k1_i
                    * k3_i
                )
                * dE3_i_dx3_i
                + (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk3_i
                    + (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * k2_i
                    * k8_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk8_i
                    * h11_i
                )
                * dE3_i_dx4_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx1_i_dk4_i
                    + (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_i
                    * k7_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_i_dk7_i
                    * h5_i
                )
                * dE3_i_dx1_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx2_i_dk4_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_i_dk7_i
                    * h5_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_i
                    * k4_i
                    * k6_i
                )
                * dE3_i_dx2_i
            )
            * k4pp_i
            - (
                (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk3_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk8_i
                    * h11_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_i
                    * k1_i
                    * k3_i
                )
                * dE2_i_dx3_i
                + (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk3_i
                    + (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * k2_i
                    * k8_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk8_i
                    * h11_i
                )
                * dE2_i_dx4_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx1_i_dk4_i
                    + (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_i
                    * k7_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_i_dk7_i
                    * h5_i
                )
                * dE2_i_dx1_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx2_i_dk4_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_i_dk7_i
                    * h5_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_i
                    * k4_i
                    * k6_i
                )
                * dE2_i_dx2_i
            )
            * k3pp_i
            - 3
            * (
                (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk3_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk8_i
                    * h11_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_i
                    * k1_i
                    * k3_i
                )
                * dE1_i_dx3_i
                + (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk3_i
                    + (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * k2_i
                    * k8_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk8_i
                    * h11_i
                )
                * dE1_i_dx4_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx1_i_dk4_i
                    + (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_i
                    * k7_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_i_dk7_i
                    * h5_i
                )
                * dE1_i_dx1_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx2_i_dk4_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_i_dk7_i
                    * h5_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_i
                    * k4_i
                    * k6_i
                )
                * dE1_i_dx2_i
            )
            * k8_i
            + 3
            * (
                (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk3_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_i_dk8_i
                    * h11_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_i
                    * k1_i
                    * k3_i
                )
                * dE4_i_dx3_i
                + (
                    (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk3_i
                    + (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * k2_i
                    * k8_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_i_dk8_i
                    * h11_i
                )
                * dE4_i_dx4_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx1_i_dk4_i
                    + (
                        wnaca
                        * (
                            dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_i_dhna * dh9_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_i
                    * k7_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_i_dk7_i
                    * h5_i
                )
                * dE4_i_dx1_i
                + (
                    (
                        wnaca
                        * (
                            dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_i_dhca * dvfrt_dv
                        + dh1_i_dhna
                        * dh3_i_dh1_i
                        * dhna_dvfrt
                        * dk4p_i_dh3_i
                        * dvfrt_dv
                    )
                    * dx2_i_dk4_i
                    + wna
                    * (
                        dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_i_dk7_i
                    * h5_i
                    + wna
                    * (
                        dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_i
                    * k4_i
                    * k6_i
                )
                * dE4_i_dx2_i
            )
            * k7_i
            + wnaca
            * (
                dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_i_dk4pp_i
            + wnaca
            * (
                dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_i_dk3pp_i
            + wna
            * (
                dh2_i_dhna * dhna_dvfrt * dvfrt_dv
                + dh1_i_dhna * dh2_i_dh1_i * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_i_dk7_i
            * h5_i
            + wna
            * (
                dh8_i_dhna * dhna_dvfrt * dvfrt_dv
                + dh7_i_dhna * dh8_i_dh7_i * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_i_dk8_i
            * h11_i
        )
        * dINaCa_i_dJncxNa_i
        - (
            (
                (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk3_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk8_ss
                    * h11_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_ss
                    * k1_ss
                    * k3_ss
                )
                * dE3_ss_dx3_ss
                + (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk3_ss
                    + (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * k2_ss
                    * k8_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk8_ss
                    * h11_ss
                )
                * dE3_ss_dx4_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_ss
                    * k7_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_ss_dk7_ss
                    * h5_ss
                )
                * dE3_ss_dx1_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx2_ss_dk4_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_ss_dk7_ss
                    * h5_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_ss
                    * k4_ss
                    * k6_ss
                )
                * dE3_ss_dx2_ss
            )
            * k4pp_ss
            - (
                (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk3_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk8_ss
                    * h11_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_ss
                    * k1_ss
                    * k3_ss
                )
                * dE2_ss_dx3_ss
                + (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk3_ss
                    + (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * k2_ss
                    * k8_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk8_ss
                    * h11_ss
                )
                * dE2_ss_dx4_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_ss
                    * k7_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_ss_dk7_ss
                    * h5_ss
                )
                * dE2_ss_dx1_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx2_ss_dk4_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_ss_dk7_ss
                    * h5_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_ss
                    * k4_ss
                    * k6_ss
                )
                * dE2_ss_dx2_ss
            )
            * k3pp_ss
            - 3
            * (
                (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk3_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk8_ss
                    * h11_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_ss
                    * k1_ss
                    * k3_ss
                )
                * dE1_ss_dx3_ss
                + (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk3_ss
                    + (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * k2_ss
                    * k8_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk8_ss
                    * h11_ss
                )
                * dE1_ss_dx4_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_ss
                    * k7_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_ss_dk7_ss
                    * h5_ss
                )
                * dE1_ss_dx1_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx2_ss_dk4_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_ss_dk7_ss
                    * h5_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_ss
                    * k4_ss
                    * k6_ss
                )
                * dE1_ss_dx2_ss
            )
            * k8_ss
            + 3
            * (
                (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk3_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx3_ss_dk8_ss
                    * h11_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h5_ss
                    * k1_ss
                    * k3_ss
                )
                * dE4_ss_dx3_ss
                + (
                    (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk3_ss
                    + (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * k2_ss
                    * k8_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx4_ss_dk8_ss
                    * h11_ss
                )
                * dE4_ss_dx4_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wnaca
                        * (
                            dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + wca * dh7_ss_dhna * dh9_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * k5_ss
                    * k7_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx1_ss_dk7_ss
                    * h5_ss
                )
                * dE4_ss_dx1_ss
                + (
                    (
                        wnaca
                        * (
                            dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                            + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                        )
                        + dhca_dvfrt * dk4p_ss_dhca * dvfrt_dv
                        + dh1_ss_dhna
                        * dh3_ss_dh1_ss
                        * dhna_dvfrt
                        * dk4p_ss_dh3_ss
                        * dvfrt_dv
                    )
                    * dx2_ss_dk4_ss
                    + wna
                    * (
                        dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * dx2_ss_dk7_ss
                    * h5_ss
                    + wna
                    * (
                        dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                        + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
                    )
                    * h11_ss
                    * k4_ss
                    * k6_ss
                )
                * dE4_ss_dx2_ss
            )
            * k7_ss
            + wnaca
            * (
                dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_ss_dk4pp_ss
            + wnaca
            * (
                dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_ss_dk3pp_ss
            + wna
            * (
                dh2_ss_dhna * dhna_dvfrt * dvfrt_dv
                + dh1_ss_dhna * dh2_ss_dh1_ss * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_ss_dk7_ss
            * h5_ss
            + wna
            * (
                dh8_ss_dhna * dhna_dvfrt * dvfrt_dv
                + dh7_ss_dhna * dh8_ss_dh7_ss * dhna_dvfrt * dvfrt_dv
            )
            * dJncxNa_ss_dk8_ss
            * h11_ss
        )
        * dINaCa_ss_dJncxNa_ss
        - dICab_dvffrt * dvffrt_dv
        - dICab_dvfrt * dvfrt_dv
        - dIKb_dxkb * dxkb_dv
        - dINab_dvffrt * dvffrt_dv
        - dINab_dvfrt * dvfrt_dv
        - zk
        * (
            -2
            * (
                (
                    dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx3_da3
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx3_db2
                )
                * dE3_dx3
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx1_da1
                    + dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx1_db4
                    + b3 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE3_dx1
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx4_da1
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx4_db2
                    + a1 * a4 * dKnao_dvfrt * da3_dKnao * dvfrt_dv
                    + b2 * b3 * dKnai_dvfrt * db4_dKnai * dvfrt_dv
                )
                * dE3_dx4
                + (
                    dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx2_db4
                    + dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx2_da3
                    + a2 * a3 * dKnai_dvfrt * da1_dKnai * dvfrt_dv
                    + b1 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE3_dx2
            )
            * a1
            + 2
            * (
                (
                    dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx3_da3
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx3_db2
                )
                * dE4_dx3
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx1_da1
                    + dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx1_db4
                    + b3 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE4_dx1
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx4_da1
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx4_db2
                    + a1 * a4 * dKnao_dvfrt * da3_dKnao * dvfrt_dv
                    + b2 * b3 * dKnai_dvfrt * db4_dKnai * dvfrt_dv
                )
                * dE4_dx4
                + (
                    dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx2_db4
                    + dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx2_da3
                    + a2 * a3 * dKnai_dvfrt * da1_dKnai * dvfrt_dv
                    + b1 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE4_dx2
            )
            * b1
            + dJnakK_da1 * dKnai_dvfrt * da1_dKnai * dvfrt_dv
        )
        * Pnak
        - zna
        * (
            -3
            * (
                (
                    dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx3_da3
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx3_db2
                )
                * dE2_dx3
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx1_da1
                    + dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx1_db4
                    + b3 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE2_dx1
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx4_da1
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx4_db2
                    + a1 * a4 * dKnao_dvfrt * da3_dKnao * dvfrt_dv
                    + b2 * b3 * dKnai_dvfrt * db4_dKnai * dvfrt_dv
                )
                * dE2_dx4
                + (
                    dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx2_db4
                    + dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx2_da3
                    + a2 * a3 * dKnai_dvfrt * da1_dKnai * dvfrt_dv
                    + b1 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE2_dx2
            )
            * b3
            + 3
            * (
                (
                    dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx3_da3
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx3_db2
                )
                * dE1_dx3
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx1_da1
                    + dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx1_db4
                    + b3 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE1_dx1
                + (
                    dKnai_dvfrt * da1_dKnai * dvfrt_dv * dx4_da1
                    + dKnao_dvfrt * db2_dKnao * dvfrt_dv * dx4_db2
                    + a1 * a4 * dKnao_dvfrt * da3_dKnao * dvfrt_dv
                    + b2 * b3 * dKnai_dvfrt * db4_dKnai * dvfrt_dv
                )
                * dE1_dx4
                + (
                    dKnai_dvfrt * db4_dKnai * dvfrt_dv * dx2_db4
                    + dKnao_dvfrt * da3_dKnao * dvfrt_dv * dx2_da3
                    + a2 * a3 * dKnai_dvfrt * da1_dKnai * dvfrt_dv
                    + b1 * b4 * dKnao_dvfrt * db2_dKnao * dvfrt_dv
                )
                * dE1_dx2
            )
            * a3
            + dJnakNa_da3 * dKnao_dvfrt * da3_dKnao * dvfrt_dv
        )
        * Pnak
    )
    states[34] = (
        np.where(
            np.abs(dv_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dv_dt_linearized)) * dv_dt / dv_dt_linearized,
            dt * dv_dt,
        )
        + v
    )

    # Expressions for the Intracellular ions component
    cmdnmax = np.where(celltype == 1, 1.3 * cmdnmax_b, cmdnmax_b)
    dnai_dt = JdiffNa * vss / vmyo + (
        -ICaNa_i - INa - INaL - INab - 3 * INaCa_i - 3 * INaK
    ) * Acap / (F * vmyo)
    dENa_dnai = -R * T / (F * zna * nai)
    dINa_dENa = -GNa * (m * m * m) * ((1 - fINap) * h * j + fINap * hp * jp)
    dINaL_dENa = -((1 - fINaLp) * hL + fINaLp * hLp) * GNaL * mL
    dINab_dnai = PNab * np.exp(vfrt) * vffrt / (-1 + np.exp(vfrt))
    dJdiffNa_dnai = -1 / tauNa
    dx1_db3 = a1 * a2 + a2 * b4 + b2 * b4
    dx4_db3 = a1 * b2 + b2 * b4
    dJnakNa_db3 = (
        -3 * E2
        - 3 * (dE2_dx1 * dx1_db3 + dE2_dx4 * dx4_db3 + b1 * b2 * dE2_dx3) * b3
        + 3 * (dE1_dx1 * dx1_db3 + dE1_dx4 * dx4_db3 + b1 * b2 * dE1_dx3) * a3
    )
    dP_dnai = -eP / (
        Knap
        * (
            (1 + H / Khp + nai / Knap + ki / Kxkur)
            * (1 + H / Khp + nai / Knap + ki / Kxkur)
        )
    )
    dPhiCaNa_i_dgamma_nai = np.exp(vfrt) * nai * vffrt / (-1 + np.exp(vfrt))
    dgamma_nai_dIi = (
        -(
            -0.3
            - 1 / (2 * ((1 + np.sqrt(Ii)) * (1 + np.sqrt(Ii))))
            + 1 / (2 * (1 + np.sqrt(Ii)) * np.sqrt(Ii))
        )
        * constA
        * np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    )
    dPhiCaNa_i_dnai = (
        (np.exp(vfrt) * gamma_nai + 0.0005 * dgamma_nai_dIi * np.exp(vfrt) * nai)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    da1_dnai = 3 * k1p * (nai * nai) / (
        (
            -1
            + ((1 + ki / Kki) * (1 + ki / Kki))
            + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
        )
        * (Knai * Knai * Knai)
    ) - 3 * k1p * ((1 + nai / Knai) * (1 + nai / Knai)) * (nai * nai * nai) / (
        (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
        )
        * np.power(Knai, 4)
    )
    db3_dP = H * k3m / (1 + MgATP / Kmgatp)
    db4_dnai = (
        -3
        * k4m
        * ((1 + nai / Knai) * (1 + nai / Knai))
        * (ki * ki)
        / (
            (Kki * Kki)
            * (
                (
                    -1
                    + ((1 + ki / Kki) * (1 + ki / Kki))
                    + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
                )
                * (
                    -1
                    + ((1 + ki / Kki) * (1 + ki / Kki))
                    + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
                )
            )
            * Knai
        )
    )
    dh1_i_dnai = (1 + hna) / kna3
    dh2_i_dnai = hna / (kna3 * h1_i) - dh1_i_dnai * hna * nai / (kna3 * (h1_i * h1_i))
    dh4_i_dnai = (1 + nai / kna2) / kna1 + nai / (kna1 * kna2)
    dh5_i_dh4_i = -(nai * nai) / (kna1 * kna2 * (h4_i * h4_i))
    dh5_i_dnai = 2 * nai / (kna1 * kna2 * h4_i) - (nai * nai) * dh4_i_dnai / (
        kna1 * kna2 * (h4_i * h4_i)
    )
    dh6_i_dh4_i = -1 / (h4_i * h4_i)
    dx2_i_dk6_i = (k1_i + k8_i) * k4_i
    dx3_i_dk6_i = (k2_i + k3_i) * k8_i + k1_i * k3_i
    dnai_dt_linearized = dJdiffNa_dnai * vss / vmyo + (
        -dINab_dnai
        - (0.0005 * dPhiCaNa_i_dgamma_nai * dgamma_nai_dIi + dPhiCaNa_i_dnai)
        * dICaNa_i_dPhiCaNa_i
        - dENa_dnai * dINaL_dENa
        - dENa_dnai * dINa_dENa
        - 3
        * (
            (
                (
                    (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * k1_i
                    * k3_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx3_i_dk6_i
                )
                * dE2_i_dx3_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx1_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx1_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * k2_i * k4_i
                )
                * dE2_i_dx1_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx2_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx2_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx2_i_dk6_i
                )
                * dE2_i_dx2_i
                + (
                    wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                    + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                )
                * dE2_i_dx4_i
                * k2_i
                * k8_i
            )
            * k2_i
            - (
                (
                    (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * k1_i
                    * k3_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx3_i_dk6_i
                )
                * dE1_i_dx3_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx1_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx1_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * k2_i * k4_i
                )
                * dE1_i_dx1_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx2_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx2_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx2_i_dk6_i
                )
                * dE1_i_dx2_i
                + (
                    wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                    + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                )
                * dE1_i_dx4_i
                * k2_i
                * k8_i
            )
            * k1_i
        )
        * dINaCa_i_dJncxCa_i
        - 3
        * (
            (
                wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
            )
            * dJncxNa_i_dk7_i
            + (
                (
                    (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * k1_i
                    * k3_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx3_i_dk6_i
                )
                * dE3_i_dx3_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx1_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx1_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * k2_i * k4_i
                )
                * dE3_i_dx1_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx2_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx2_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx2_i_dk6_i
                )
                * dE3_i_dx2_i
                + (
                    wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                    + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                )
                * dE3_i_dx4_i
                * k2_i
                * k8_i
            )
            * k4pp_i
            - (
                (
                    (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * k1_i
                    * k3_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx3_i_dk6_i
                )
                * dE2_i_dx3_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx1_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx1_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * k2_i * k4_i
                )
                * dE2_i_dx1_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx2_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx2_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx2_i_dk6_i
                )
                * dE2_i_dx2_i
                + (
                    wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                    + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                )
                * dE2_i_dx4_i
                * k2_i
                * k8_i
            )
            * k3pp_i
            - 3
            * (
                (
                    (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * k1_i
                    * k3_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx3_i_dk6_i
                )
                * dE1_i_dx3_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx1_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx1_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * k2_i * k4_i
                )
                * dE1_i_dx1_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx2_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx2_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx2_i_dk6_i
                )
                * dE1_i_dx2_i
                + (
                    wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                    + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                )
                * dE1_i_dx4_i
                * k2_i
                * k8_i
            )
            * k8_i
            + 3
            * (
                (
                    (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * k1_i
                    * k3_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx3_i_dk6_i
                )
                * dE4_i_dx3_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx1_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx1_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * k2_i * k4_i
                )
                * dE4_i_dx1_i
                + (
                    (
                        wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                        + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                    )
                    * dx2_i_dk4_i
                    + (
                        wna * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * h5_i
                        + wna * (dh4_i_dnai * dh5_i_dh4_i + dh5_i_dnai) * h2_i
                    )
                    * dx2_i_dk7_i
                    + kcaon * cai * dh4_i_dnai * dh6_i_dh4_i * dx2_i_dk6_i
                )
                * dE4_i_dx2_i
                + (
                    wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai)
                    + dh1_i_dnai * dh3_i_dh1_i * dk4p_i_dh3_i
                )
                * dE4_i_dx4_i
                * k2_i
                * k8_i
            )
            * k7_i
            + wnaca * (dh1_i_dnai * dh2_i_dh1_i + dh2_i_dnai) * dJncxNa_i_dk4pp_i
        )
        * dINaCa_i_dJncxNa_i
        - 3
        * zk
        * (
            dJnakK_da1 * da1_dnai
            - 2
            * (
                (db4_dnai * dx2_db4 + a2 * a3 * da1_dnai) * dE3_dx2
                + (da1_dnai * dx1_da1 + db4_dnai * dx1_db4 + dP_dnai * db3_dP * dx1_db3)
                * dE3_dx1
                + (da1_dnai * dx4_da1 + b2 * b3 * db4_dnai + dP_dnai * db3_dP * dx4_db3)
                * dE3_dx4
                + b1 * b2 * dE3_dx3 * dP_dnai * db3_dP
            )
            * a1
            + 2
            * (
                (db4_dnai * dx2_db4 + a2 * a3 * da1_dnai) * dE4_dx2
                + (da1_dnai * dx1_da1 + db4_dnai * dx1_db4 + dP_dnai * db3_dP * dx1_db3)
                * dE4_dx1
                + (da1_dnai * dx4_da1 + b2 * b3 * db4_dnai + dP_dnai * db3_dP * dx4_db3)
                * dE4_dx4
                + b1 * b2 * dE4_dx3 * dP_dnai * db3_dP
            )
            * b1
        )
        * Pnak
        - 3
        * zna
        * (
            -3
            * (
                (db4_dnai * dx2_db4 + a2 * a3 * da1_dnai) * dE2_dx2
                + (da1_dnai * dx1_da1 + db4_dnai * dx1_db4 + dP_dnai * db3_dP * dx1_db3)
                * dE2_dx1
                + (da1_dnai * dx4_da1 + b2 * b3 * db4_dnai + dP_dnai * db3_dP * dx4_db3)
                * dE2_dx4
                + b1 * b2 * dE2_dx3 * dP_dnai * db3_dP
            )
            * b3
            + 3
            * (
                (db4_dnai * dx2_db4 + a2 * a3 * da1_dnai) * dE1_dx2
                + (da1_dnai * dx1_da1 + db4_dnai * dx1_db4 + dP_dnai * db3_dP * dx1_db3)
                * dE1_dx1
                + (da1_dnai * dx4_da1 + b2 * b3 * db4_dnai + dP_dnai * db3_dP * dx4_db3)
                * dE1_dx4
                + b1 * b2 * dE1_dx3 * dP_dnai * db3_dP
            )
            * a3
            + dJnakNa_db3 * dP_dnai * db3_dP
        )
        * Pnak
    ) * Acap / (F * vmyo)
    states[35] = (
        np.where(
            np.abs(dnai_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dnai_dt_linearized)) * dnai_dt / dnai_dt_linearized,
            dt * dnai_dt,
        )
        + nai
    )
    dnass_dt = -JdiffNa + (-ICaNa_ss - 3 * INaCa_ss) * Acap / (F * vss)
    dPhiCaNa_ss_dgamma_nass = np.exp(vfrt) * nass * vffrt / (-1 + np.exp(vfrt))
    dgamma_nass_dIss = (
        -(
            -0.3
            - 1 / (2 * ((1 + np.sqrt(Iss)) * (1 + np.sqrt(Iss))))
            + 1 / (2 * (1 + np.sqrt(Iss)) * np.sqrt(Iss))
        )
        * constA
        * np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    )
    dPhiCaNa_ss_dnass = (
        (np.exp(vfrt) * gamma_nass + 0.0005 * dgamma_nass_dIss * np.exp(vfrt) * nass)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    dh1_ss_dnass = (1 + hna) / kna3
    dh2_ss_dnass = hna / (kna3 * h1_ss) - dh1_ss_dnass * hna * nass / (
        kna3 * (h1_ss * h1_ss)
    )
    dh4_ss_dnass = (1 + nass / kna2) / kna1 + nass / (kna1 * kna2)
    dh5_ss_dh4_ss = -(nass * nass) / (kna1 * kna2 * (h4_ss * h4_ss))
    dh5_ss_dnass = 2 * nass / (kna1 * kna2 * h4_ss) - (nass * nass) * dh4_ss_dnass / (
        kna1 * kna2 * (h4_ss * h4_ss)
    )
    dh6_ss_dh4_ss = -1 / (h4_ss * h4_ss)
    dx2_ss_dk6_ss = (k1_ss + k8_ss) * k4_ss
    dx3_ss_dk6_ss = (k2_ss + k3_ss) * k8_ss + k1_ss * k3_ss
    dnass_dt_linearized = -1 / tauNa + (
        -(0.0005 * dPhiCaNa_ss_dgamma_nass * dgamma_nass_dIss + dPhiCaNa_ss_dnass)
        * dICaNa_ss_dPhiCaNa_ss
        - 3
        * (
            (
                (
                    (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * k1_ss
                    * k3_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx3_ss_dk6_ss
                )
                * dE2_ss_dx3_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx1_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * k2_ss * k4_ss
                )
                * dE2_ss_dx1_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx2_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx2_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx2_ss_dk6_ss
                )
                * dE2_ss_dx2_ss
                + (
                    wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                    + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                )
                * dE2_ss_dx4_ss
                * k2_ss
                * k8_ss
            )
            * k2_ss
            - (
                (
                    (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * k1_ss
                    * k3_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx3_ss_dk6_ss
                )
                * dE1_ss_dx3_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx1_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * k2_ss * k4_ss
                )
                * dE1_ss_dx1_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx2_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx2_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx2_ss_dk6_ss
                )
                * dE1_ss_dx2_ss
                + (
                    wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                    + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                )
                * dE1_ss_dx4_ss
                * k2_ss
                * k8_ss
            )
            * k1_ss
        )
        * dINaCa_ss_dJncxCa_ss
        - 3
        * (
            (
                wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
            )
            * dJncxNa_ss_dk7_ss
            + (
                (
                    (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * k1_ss
                    * k3_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx3_ss_dk6_ss
                )
                * dE3_ss_dx3_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx1_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * k2_ss * k4_ss
                )
                * dE3_ss_dx1_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx2_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx2_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx2_ss_dk6_ss
                )
                * dE3_ss_dx2_ss
                + (
                    wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                    + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                )
                * dE3_ss_dx4_ss
                * k2_ss
                * k8_ss
            )
            * k4pp_ss
            - (
                (
                    (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * k1_ss
                    * k3_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx3_ss_dk6_ss
                )
                * dE2_ss_dx3_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx1_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * k2_ss * k4_ss
                )
                * dE2_ss_dx1_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx2_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx2_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx2_ss_dk6_ss
                )
                * dE2_ss_dx2_ss
                + (
                    wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                    + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                )
                * dE2_ss_dx4_ss
                * k2_ss
                * k8_ss
            )
            * k3pp_ss
            - 3
            * (
                (
                    (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * k1_ss
                    * k3_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx3_ss_dk6_ss
                )
                * dE1_ss_dx3_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx1_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * k2_ss * k4_ss
                )
                * dE1_ss_dx1_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx2_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx2_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx2_ss_dk6_ss
                )
                * dE1_ss_dx2_ss
                + (
                    wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                    + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                )
                * dE1_ss_dx4_ss
                * k2_ss
                * k8_ss
            )
            * k8_ss
            + 3
            * (
                (
                    (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * k1_ss
                    * k3_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx3_ss_dk6_ss
                )
                * dE4_ss_dx3_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx1_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx1_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * k2_ss * k4_ss
                )
                * dE4_ss_dx1_ss
                + (
                    (
                        wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                        + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                    )
                    * dx2_ss_dk4_ss
                    + (
                        wna * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass) * h5_ss
                        + wna * (dh4_ss_dnass * dh5_ss_dh4_ss + dh5_ss_dnass) * h2_ss
                    )
                    * dx2_ss_dk7_ss
                    + kcaon * cass * dh4_ss_dnass * dh6_ss_dh4_ss * dx2_ss_dk6_ss
                )
                * dE4_ss_dx2_ss
                + (
                    wnaca * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
                    + dh1_ss_dnass * dh3_ss_dh1_ss * dk4p_ss_dh3_ss
                )
                * dE4_ss_dx4_ss
                * k2_ss
                * k8_ss
            )
            * k7_ss
            + wnaca
            * (dh1_ss_dnass * dh2_ss_dh1_ss + dh2_ss_dnass)
            * dJncxNa_ss_dk4pp_ss
        )
        * dINaCa_ss_dJncxNa_ss
    ) * Acap / (F * vss)
    states[36] = (
        np.where(
            np.abs(dnass_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dnass_dt_linearized)) * dnass_dt / dnass_dt_linearized,
            dt * dnass_dt,
        )
        + nass
    )
    dki_dt = JdiffK * vss / vmyo + (
        -ICaK_i - IK1 - IKb - IKr - IKs - I_katp - Istim - Ito + 2 * INaK
    ) * Acap / (F * vmyo)
    dEK_dki = -R * T / (F * zk * ki)
    dEKs_dki = -R * T / (F * zk * (PKNa * nai + ki))
    daK1_dEK = (
        0.0011435228161993972
        * np.exp(0.1217 * v - 0.1217 * EK)
        / (
            (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
            * (1 + 0.0022951253918281865 * np.exp(0.1217 * v - 0.1217 * EK))
        )
    )
    dbK1_dEK = (
        -6.919349479611957e-18 * np.exp(0.0618 * v - 0.0618 * EK)
        - 0.8506978434250505 * np.exp(0.0674 * v - 0.0674 * EK)
    ) / (
        1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v)
    ) - 0.016099950792318272 * (
        12.621629724407278 * np.exp(0.0674 * v - 0.0674 * EK)
        + 1.1196358381249121e-16 * np.exp(0.0618 * v - 0.0618 * EK)
    ) * np.exp(
        0.1629 * EK - 0.1629 * v
    ) / (
        (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
        * (1 + 0.09883333819716558 * np.exp(0.1629 * EK - 0.1629 * v))
    )
    dIK1_dEK = (
        -np.sqrt(5) * np.sqrt(ko) * GK1 * K1ss / 5
        + np.sqrt(5)
        * np.sqrt(ko)
        * (-EK + v)
        * (dK1ss_daK1 * daK1_dEK + dK1ss_dbK1 * dbK1_dEK)
        * GK1
        / 5
    )
    dIKb_dEK = -GKb * xkb
    dIKr_dEK = -np.sqrt(5) * np.sqrt(ko) * GKr * O / 5
    dIKs_dEKs = -GKs * KsCa * xs1 * xs2
    dI_katp_dEK = -fkatp * gkatp * akik * bkik
    dIto_dEK = -((1 - fItop) * a * i + ap * fItop * ip) * Gto
    dJdiffK_dki = -1 / tauK
    dP_dki = -eP / (
        Kxkur
        * (
            (1 + H / Khp + nai / Knap + ki / Kxkur)
            * (1 + H / Khp + nai / Knap + ki / Kxkur)
        )
    )
    dPhiCaK_i_dgamma_ki = np.exp(vfrt) * ki * vffrt / (-1 + np.exp(vfrt))
    dgamma_ki_dIi = (
        -(
            -0.3
            - 1 / (2 * ((1 + np.sqrt(Ii)) * (1 + np.sqrt(Ii))))
            + 1 / (2 * (1 + np.sqrt(Ii)) * np.sqrt(Ii))
        )
        * constA
        * np.exp(-(-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    )
    dPhiCaK_i_dki = (
        (np.exp(vfrt) * gamma_ki + 0.0005 * dgamma_ki_dIi * np.exp(vfrt) * ki)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    da1_dki = (
        -2
        * k1p
        * (nai * nai * nai)
        * (1 + ki / Kki)
        / (
            Kki
            * (
                (
                    -1
                    + ((1 + ki / Kki) * (1 + ki / Kki))
                    + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
                )
                * (
                    -1
                    + ((1 + ki / Kki) * (1 + ki / Kki))
                    + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
                )
            )
            * (Knai * Knai * Knai)
        )
    )
    db4_dki = 2 * k4m * ki / (
        (Kki * Kki)
        * (
            -1
            + ((1 + ki / Kki) * (1 + ki / Kki))
            + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
        )
    ) - 2 * k4m * (ki * ki) * (1 + ki / Kki) / (
        (Kki * Kki * Kki)
        * (
            (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
            * (
                -1
                + ((1 + ki / Kki) * (1 + ki / Kki))
                + ((1 + nai / Knai) * (1 + nai / Knai) * (1 + nai / Knai))
            )
        )
    )
    dki_dt_linearized = dJdiffK_dki * vss / vmyo + (
        -(0.0005 * dPhiCaK_i_dgamma_ki * dgamma_ki_dIi + dPhiCaK_i_dki)
        * dICaK_i_dPhiCaK_i
        - (dEK_dki * dK1ss_daK1 * daK1_dEK + dEK_dki * dK1ss_dbK1 * dbK1_dEK)
        * dIK1_dK1ss
        - dEK_dki * dIK1_dEK
        - dEK_dki * dIKb_dEK
        - dEK_dki * dIKr_dEK
        - dEK_dki * dI_katp_dEK
        - dEK_dki * dIto_dEK
        - dEKs_dki * dIKs_dEKs
        + 2
        * zk
        * (
            dJnakK_da1 * da1_dki
            - 2
            * (
                (db4_dki * dx2_db4 + a2 * a3 * da1_dki) * dE3_dx2
                + (da1_dki * dx1_da1 + db4_dki * dx1_db4 + dP_dki * db3_dP * dx1_db3)
                * dE3_dx1
                + (da1_dki * dx4_da1 + b2 * b3 * db4_dki + dP_dki * db3_dP * dx4_db3)
                * dE3_dx4
                + b1 * b2 * dE3_dx3 * dP_dki * db3_dP
            )
            * a1
            + 2
            * (
                (db4_dki * dx2_db4 + a2 * a3 * da1_dki) * dE4_dx2
                + (da1_dki * dx1_da1 + db4_dki * dx1_db4 + dP_dki * db3_dP * dx1_db3)
                * dE4_dx1
                + (da1_dki * dx4_da1 + b2 * b3 * db4_dki + dP_dki * db3_dP * dx4_db3)
                * dE4_dx4
                + b1 * b2 * dE4_dx3 * dP_dki * db3_dP
            )
            * b1
        )
        * Pnak
        + 2
        * zna
        * (
            -3
            * (
                (db4_dki * dx2_db4 + a2 * a3 * da1_dki) * dE2_dx2
                + (da1_dki * dx1_da1 + db4_dki * dx1_db4 + dP_dki * db3_dP * dx1_db3)
                * dE2_dx1
                + (da1_dki * dx4_da1 + b2 * b3 * db4_dki + dP_dki * db3_dP * dx4_db3)
                * dE2_dx4
                + b1 * b2 * dE2_dx3 * dP_dki * db3_dP
            )
            * b3
            + 3
            * (
                (db4_dki * dx2_db4 + a2 * a3 * da1_dki) * dE1_dx2
                + (da1_dki * dx1_da1 + db4_dki * dx1_db4 + dP_dki * db3_dP * dx1_db3)
                * dE1_dx1
                + (da1_dki * dx4_da1 + b2 * b3 * db4_dki + dP_dki * db3_dP * dx4_db3)
                * dE1_dx4
                + b1 * b2 * dE1_dx3 * dP_dki * db3_dP
            )
            * a3
            + dJnakNa_db3 * dP_dki * db3_dP
        )
        * Pnak
    ) * Acap / (F * vmyo)
    states[37] = (
        np.where(
            np.abs(dki_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dki_dt_linearized)) * dki_dt / dki_dt_linearized,
            dt * dki_dt,
        )
        + ki
    )
    dkss_dt = -JdiffK - Acap * ICaK_ss / (F * vss)
    dPhiCaK_ss_dgamma_kss = np.exp(vfrt) * kss * vffrt / (-1 + np.exp(vfrt))
    dgamma_kss_dIss = (
        -(
            -0.3
            - 1 / (2 * ((1 + np.sqrt(Iss)) * (1 + np.sqrt(Iss))))
            + 1 / (2 * (1 + np.sqrt(Iss)) * np.sqrt(Iss))
        )
        * constA
        * np.exp(-(-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    )
    dPhiCaK_ss_dkss = (
        (np.exp(vfrt) * gamma_kss + 0.0005 * dgamma_kss_dIss * np.exp(vfrt) * kss)
        * vffrt
        / (-1 + np.exp(vfrt))
    )
    dkss_dt_linearized = -1 / tauK - (
        0.0005 * dPhiCaK_ss_dgamma_kss * dgamma_kss_dIss + dPhiCaK_ss_dkss
    ) * Acap * dICaK_ss_dPhiCaK_ss / (F * vss)
    states[38] = (
        np.where(
            np.abs(dkss_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dkss_dt_linearized)) * dkss_dt / dkss_dt_linearized,
            dt * dkss_dt,
        )
        + kss
    )
    dcli_dt = JdiffCl * vss / vmyo + (IClCa_sl + IClb) * Acap / (F * vmyo)
    dECl_dcli = -R * T / (F * zcl * cli)
    dIClCa_sl_dECl = -GClCa * (1 - Fjunc) / (1 + KdClCa / cai)
    dJdiffCl_dcli = -1 / tauNa
    dcli_dt_linearized = dJdiffCl_dcli * vss / vmyo + (
        dECl_dcli * dIClCa_sl_dECl - GClb * dECl_dcli
    ) * Acap / (F * vmyo)
    states[39] = (
        np.where(
            np.abs(dcli_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dcli_dt_linearized)) * dcli_dt / dcli_dt_linearized,
            dt * dcli_dt,
        )
        + cli
    )
    dclss_dt = -JdiffCl + Acap * IClCa_junc / (F * vss)
    dEClss_dclss = -R * T / (F * zcl * clss)
    dIClCa_junc_dEClss = -Fjunc * GClCa / (1 + KdClCa / cass)
    dclss_dt_linearized = -1 / tauNa + Acap * dEClss_dclss * dIClCa_junc_dEClss / (
        F * vss
    )
    states[40] = (
        np.where(
            np.abs(dclss_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dclss_dt_linearized)) * dclss_dt / dclss_dt_linearized,
            dt * dclss_dt,
        )
        + clss
    )
    Bcai = 1.0 / (
        1
        + kmcmdn * cmdnmax / ((kmcmdn + cai) * (kmcmdn + cai))
        + kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai))
    )
    dcai_dt = (
        Jdiff * vss / vmyo
        - Jup * vnsr / vmyo
        + (-ICaL_i - ICab - IpCa + 2 * INaCa_i) * Acap / (2 * F * vmyo)
    ) * Bcai
    dBcai_dcai = (
        2 * kmcmdn * cmdnmax / ((kmcmdn + cai) * (kmcmdn + cai) * (kmcmdn + cai))
        + 2 * kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai) * (kmtrpn + cai))
    ) / (
        (
            1
            + kmcmdn * cmdnmax / ((kmcmdn + cai) * (kmcmdn + cai))
            + kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai))
        )
        * (
            1
            + kmcmdn * cmdnmax / ((kmcmdn + cai) * (kmcmdn + cai))
            + kmtrpn * trpnmax / ((kmtrpn + cai) * (kmtrpn + cai))
        )
    )
    dgamma_cai_dIi = (
        -4
        * (
            -0.3
            - 1 / (2 * ((1 + np.sqrt(Ii)) * (1 + np.sqrt(Ii))))
            + 1 / (2 * (1 + np.sqrt(Ii)) * np.sqrt(Ii))
        )
        * constA
        * np.exp(-4 * (-0.3 * Ii + np.sqrt(Ii) / (1 + np.sqrt(Ii))) * constA)
    )
    dICab_dcai = (
        4
        * PCab
        * (
            np.exp(2 * vfrt) * gamma_cai
            + 0.002 * cai * dgamma_cai_dIi * np.exp(2 * vfrt)
        )
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    dICab_dgamma_cai = (
        4 * PCab * cai * np.exp(2 * vfrt) * vffrt / (-1 + np.exp(2 * vfrt))
    )
    dINaCa_i_dallo_i = (1 - INaCa_fractionSS) * (zca * JncxCa_i + zna * JncxNa_i) * Gncx
    dIpCa_dcai = GpCa / (KmCap + cai) - GpCa * cai / ((KmCap + cai) * (KmCap + cai))
    dJdiff_dcai = -1 / tauCa
    dJup_dJupnp = Jup_b * (1 - fJupp)
    dJupnp_dcai = 0.005425 * upScale / (0.00092 + cai) - 0.005425 * cai * upScale / (
        (0.00092 + cai) * (0.00092 + cai)
    )
    dJupp_dcai = 0.01491875 * upScale / (0.00075 + cai) - 0.01491875 * cai * upScale / (
        (0.00075 + cai) * (0.00075 + cai)
    )
    dPhiCaL_i_dcai = (
        4
        * (
            np.exp(2 * vfrt) * gamma_cai
            + 0.002 * cai * dgamma_cai_dIi * np.exp(2 * vfrt)
        )
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    dPhiCaL_i_dgamma_cai = 4 * cai * np.exp(2 * vfrt) * vffrt / (-1 + np.exp(2 * vfrt))
    dallo_i_dcai = (
        2
        * (KmCaAct * KmCaAct)
        / (
            (
                (1 + (KmCaAct * KmCaAct) / (cai * cai))
                * (1 + (KmCaAct * KmCaAct) / (cai * cai))
            )
            * (cai * cai * cai)
        )
    )
    dcai_dt_linearized = (
        Jdiff * vss / vmyo
        - Jup * vnsr / vmyo
        + (-ICaL_i - ICab - IpCa + 2 * INaCa_i) * Acap / (2 * F * vmyo)
    ) * dBcai_dcai + (
        dJdiff_dcai * vss / vmyo
        - (dJup_dJupnp * dJupnp_dcai + Jup_b * dJupp_dcai * fJupp) * vnsr / vmyo
        + (
            -dICab_dcai
            - dIpCa_dcai
            - (0.002 * dPhiCaL_i_dgamma_cai * dgamma_cai_dIi + dPhiCaL_i_dcai)
            * dICaL_i_dPhiCaL_i
            + 2
            * (
                (
                    kcaon * dE2_i_dx2_i * dx2_i_dk6_i * h6_i
                    + kcaon * dE2_i_dx3_i * dx3_i_dk6_i * h6_i
                    + kcaon * dE2_i_dx1_i * h6_i * k2_i * k4_i
                )
                * k2_i
                - (
                    kcaon * dE1_i_dx2_i * dx2_i_dk6_i * h6_i
                    + kcaon * dE1_i_dx3_i * dx3_i_dk6_i * h6_i
                    + kcaon * dE1_i_dx1_i * h6_i * k2_i * k4_i
                )
                * k1_i
            )
            * dINaCa_i_dJncxCa_i
            + 2
            * (
                (
                    kcaon * dE3_i_dx2_i * dx2_i_dk6_i * h6_i
                    + kcaon * dE3_i_dx3_i * dx3_i_dk6_i * h6_i
                    + kcaon * dE3_i_dx1_i * h6_i * k2_i * k4_i
                )
                * k4pp_i
                - (
                    kcaon * dE2_i_dx2_i * dx2_i_dk6_i * h6_i
                    + kcaon * dE2_i_dx3_i * dx3_i_dk6_i * h6_i
                    + kcaon * dE2_i_dx1_i * h6_i * k2_i * k4_i
                )
                * k3pp_i
                - 3
                * (
                    kcaon * dE1_i_dx2_i * dx2_i_dk6_i * h6_i
                    + kcaon * dE1_i_dx3_i * dx3_i_dk6_i * h6_i
                    + kcaon * dE1_i_dx1_i * h6_i * k2_i * k4_i
                )
                * k8_i
                + 3
                * (
                    kcaon * dE4_i_dx2_i * dx2_i_dk6_i * h6_i
                    + kcaon * dE4_i_dx3_i * dx3_i_dk6_i * h6_i
                    + kcaon * dE4_i_dx1_i * h6_i * k2_i * k4_i
                )
                * k7_i
            )
            * dINaCa_i_dJncxNa_i
            + 2 * dINaCa_i_dallo_i * dallo_i_dcai
            - 0.002 * dICab_dgamma_cai * dgamma_cai_dIi
        )
        * Acap
        / (2 * F * vmyo)
    ) * Bcai
    states[41] = (
        np.where(
            np.abs(dcai_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dcai_dt_linearized)) * dcai_dt / dcai_dt_linearized,
            dt * dcai_dt,
        )
        + cai
    )
    Bcass = 1.0 / (
        1
        + BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass))
        + BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass))
    )
    dcass_dt = (
        -Jdiff + Jrel * vjsr / vss + (-ICaL_ss + 2 * INaCa_ss) * Acap / (2 * F * vss)
    ) * Bcass
    dBcass_dcass = (
        2 * BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass) * (KmBSL + cass))
        + 2 * BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass) * (KmBSR + cass))
    ) / (
        (
            1
            + BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass))
            + BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass))
        )
        * (
            1
            + BSLmax * KmBSL / ((KmBSL + cass) * (KmBSL + cass))
            + BSRmax * KmBSR / ((KmBSR + cass) * (KmBSR + cass))
        )
    )
    dCaMKb_dcass = (
        CaMKo
        * KmCaM
        * (1 - CaMKt)
        / (((1 + KmCaM / cass) * (1 + KmCaM / cass)) * (cass * cass))
    )
    dICaL_ss_dfICaLp = ICaL_fractionSS * (
        ((1 - nca_ss) * fp + fcap * jca * nca_ss) * PCap * PhiCaL_ss * d
        - ((1 - nca_ss) * f + fca * jca * nca_ss) * PCa * PhiCaL_ss * d
    )
    dINaCa_ss_dallo_ss = INaCa_fractionSS * (zca * JncxCa_ss + zna * JncxNa_ss) * Gncx
    dJrel_dfJrelp = Jrel_b * (-Jrel_np + Jrel_p)
    dgamma_cass_dIss = (
        -4
        * (
            -0.3
            - 1 / (2 * ((1 + np.sqrt(Iss)) * (1 + np.sqrt(Iss))))
            + 1 / (2 * (1 + np.sqrt(Iss)) * np.sqrt(Iss))
        )
        * constA
        * np.exp(-4 * (-0.3 * Iss + np.sqrt(Iss) / (1 + np.sqrt(Iss))) * constA)
    )
    dPhiCaL_ss_dcass = (
        4
        * (
            np.exp(2 * vfrt) * gamma_cass
            + 0.002 * cass * dgamma_cass_dIss * np.exp(2 * vfrt)
        )
        * vffrt
        / (-1 + np.exp(2 * vfrt))
    )
    dPhiCaL_ss_dgamma_cass = (
        4 * cass * np.exp(2 * vfrt) * vffrt / (-1 + np.exp(2 * vfrt))
    )
    dallo_ss_dcass = (
        2
        * (KmCaAct * KmCaAct)
        / (
            (
                (1 + (KmCaAct * KmCaAct) / (cass * cass))
                * (1 + (KmCaAct * KmCaAct) / (cass * cass))
            )
            * (cass * cass * cass)
        )
    )
    dfICaLp_dCaMKa = KmCaMK / (
        ((1 + KmCaMK / CaMKa) * (1 + KmCaMK / CaMKa)) * (CaMKa * CaMKa)
    )
    dfJrelp_dCaMKa = KmCaMK / (
        ((1 + KmCaMK / CaMKa) * (1 + KmCaMK / CaMKa)) * (CaMKa * CaMKa)
    )
    dcass_dt_linearized = (
        -1 / tauCa
        + (
            -(0.002 * dPhiCaL_ss_dgamma_cass * dgamma_cass_dIss + dPhiCaL_ss_dcass)
            * dICaL_ss_dPhiCaL_ss
            + 2
            * (
                (
                    kcaon * dE2_ss_dx2_ss * dx2_ss_dk6_ss * h6_ss
                    + kcaon * dE2_ss_dx3_ss * dx3_ss_dk6_ss * h6_ss
                    + kcaon * dE2_ss_dx1_ss * h6_ss * k2_ss * k4_ss
                )
                * k2_ss
                - (
                    kcaon * dE1_ss_dx2_ss * dx2_ss_dk6_ss * h6_ss
                    + kcaon * dE1_ss_dx3_ss * dx3_ss_dk6_ss * h6_ss
                    + kcaon * dE1_ss_dx1_ss * h6_ss * k2_ss * k4_ss
                )
                * k1_ss
            )
            * dINaCa_ss_dJncxCa_ss
            + 2
            * (
                (
                    kcaon * dE3_ss_dx2_ss * dx2_ss_dk6_ss * h6_ss
                    + kcaon * dE3_ss_dx3_ss * dx3_ss_dk6_ss * h6_ss
                    + kcaon * dE3_ss_dx1_ss * h6_ss * k2_ss * k4_ss
                )
                * k4pp_ss
                - (
                    kcaon * dE2_ss_dx2_ss * dx2_ss_dk6_ss * h6_ss
                    + kcaon * dE2_ss_dx3_ss * dx3_ss_dk6_ss * h6_ss
                    + kcaon * dE2_ss_dx1_ss * h6_ss * k2_ss * k4_ss
                )
                * k3pp_ss
                - 3
                * (
                    kcaon * dE1_ss_dx2_ss * dx2_ss_dk6_ss * h6_ss
                    + kcaon * dE1_ss_dx3_ss * dx3_ss_dk6_ss * h6_ss
                    + kcaon * dE1_ss_dx1_ss * h6_ss * k2_ss * k4_ss
                )
                * k8_ss
                + 3
                * (
                    kcaon * dE4_ss_dx2_ss * dx2_ss_dk6_ss * h6_ss
                    + kcaon * dE4_ss_dx3_ss * dx3_ss_dk6_ss * h6_ss
                    + kcaon * dE4_ss_dx1_ss * h6_ss * k2_ss * k4_ss
                )
                * k7_ss
            )
            * dINaCa_ss_dJncxNa_ss
            + 2 * dINaCa_ss_dallo_ss * dallo_ss_dcass
            - dCaMKb_dcass * dICaL_ss_dfICaLp * dfICaLp_dCaMKa
        )
        * Acap
        / (2 * F * vss)
        + dCaMKb_dcass * dJrel_dfJrelp * dfJrelp_dCaMKa * vjsr / vss
    ) * Bcass + (
        -Jdiff + Jrel * vjsr / vss + (-ICaL_ss + 2 * INaCa_ss) * Acap / (2 * F * vss)
    ) * dBcass_dcass
    states[42] = (
        np.where(
            np.abs(dcass_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dcass_dt_linearized)) * dcass_dt / dcass_dt_linearized,
            dt * dcass_dt,
        )
        + cass
    )
    dcansr_dt = -Jtr * vjsr / vnsr + Jup
    dcansr_dt_linearized = -0.0003255 * Jup_b - vjsr / (60 * vnsr)
    states[43] = (
        np.where(
            np.abs(dcansr_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dcansr_dt_linearized)) * dcansr_dt / dcansr_dt_linearized,
            dt * dcansr_dt,
        )
        + cansr
    )
    Bcajsr = 1.0 / (1 + csqnmax * kmcsqn / ((kmcsqn + cajsr) * (kmcsqn + cajsr)))
    dcajsr_dt = (-Jrel + Jtr) * Bcajsr
    dBcajsr_dcajsr = (
        2
        * csqnmax
        * kmcsqn
        / (
            (
                (1 + csqnmax * kmcsqn / ((kmcsqn + cajsr) * (kmcsqn + cajsr)))
                * (1 + csqnmax * kmcsqn / ((kmcsqn + cajsr) * (kmcsqn + cajsr)))
            )
            * ((kmcsqn + cajsr) * (kmcsqn + cajsr) * (kmcsqn + cajsr))
        )
    )
    dcajsr_dt_linearized = -Bcajsr / 60 + (-Jrel + Jtr) * dBcajsr_dcajsr
    states[44] = (
        np.where(
            np.abs(dcajsr_dt_linearized) > 1e-08,
            (-1 + np.exp(dt * dcajsr_dt_linearized)) * dcajsr_dt / dcajsr_dt_linearized,
            dt * dcajsr_dt,
        )
        + cajsr
    )

    # Return results
    return states
