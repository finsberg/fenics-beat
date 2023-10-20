import dolfin
import beat
import numpy as np
import pytest


def test_ecg_recovery():
    mesh = dolfin.UnitCubeMesh(5, 5, 5)
    cfun = dolfin.MeshFunction("size_t", mesh, 3)

    class Heart(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return 0.25 < x[0] < 0.75 and 0.25 < x[1] < 0.75 and 0.25 < x[2] < 0.75

    heart = Heart()
    cfun.set_all(0)
    heart.mark(cfun, 1)
    dx = dolfin.dx(subdomain_data=cfun, domain=mesh)(1)

    # Need second order here, otherwise it will be zero
    V = dolfin.FunctionSpace(mesh, "CG", 2)
    v = dolfin.Function(V)
    v.interpolate(dolfin.Expression("100 * x[0] * x[0]", degree=2))

    u_e = beat.ecg.ecg_recovery(
        v=v,
        mesh=mesh,
        point=[1.0, 1.0, 1.0],
        dx=dx,
        sigma_b=1.0,
    )

    assert np.isclose(float(u_e), -0.04901946833095499)


def test_12_leads_ecg_no_precordial():
    N = 10
    x = np.ones(N)
    la = 1.2
    ra = 4.5
    ll = 3.6

    ecg = beat.ecg.Leads12(LA=la * x, RA=ra * x, LL=ll * x)
    assert np.allclose(ecg.Vw, np.mean([la, ra, ll]))
    assert np.allclose(ecg.aVR, ra - 0.5 * (la + ll))
    assert np.allclose(ecg.aVL, la - 0.5 * (ra + ll))
    assert np.allclose(ecg.aVF, ll - 0.5 * (la + ra))

    with pytest.raises(AttributeError):
        ecg.V1_


def test_12_leads_ecg():
    N = 10
    x = np.ones(N)
    la = 1.2
    ra = 4.5
    ll = 3.6
    v1 = 1.0
    v2 = 2.0
    v3 = 3.0
    v4 = 4.0
    v5 = 5.0
    v6 = 6.0

    Vw = np.mean([la, ra, ll])

    ecg = beat.ecg.Leads12(
        LA=la * x,
        RA=ra * x,
        LL=ll * x,
        V1=v1 * x,
        V2=v2 * x,
        V3=v3 * x,
        V4=v4 * x,
        V5=v5 * x,
        V6=v6 * x,
    )

    for i, vi in enumerate([v1, v2, v3, v4, v5, v6], start=1):
        assert np.allclose(getattr(ecg, f"V{i}_"), vi - Vw)
