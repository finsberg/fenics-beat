import dolfin
import beat
import numpy as np


def test_single_stimulation():
    mesh = dolfin.UnitIntervalMesh(10)
    value = 2.0
    end = 1.0
    start = 0.5
    dt = 0.01
    time = dolfin.Constant(0.0)
    expr = dolfin.Expression(
        "(t >= start && t <= end )? value : 0.0",
        value=value,
        start=start,
        end=end,
        t=time,
        degree=0,
    )
    I_s = beat.stimulation.Stimulus(dz=dolfin.dx, expr=expr)

    pde = beat.MonodomainModel(
        time=time,
        mesh=mesh,
        M=dolfin.Constant(0.0),
        I_s=I_s,
        params=dict(theta=1.0),
    )

    pde.step((0.0, 0.4))
    assert np.allclose(pde.state.vector().get_local(), 0.0)

    t0 = 0.9
    stim_t0 = value * (t0 - start)
    pde.solve((0.4, t0), dt=dt)
    # At time dt the stimulus should be value and since M is zero the state should be value * dt
    assert np.allclose(pde.state.vector().get_local(), stim_t0)

    pde.solve((t0, end + dt), dt=dt)

    # At end the stimulus should be zero and since M is zero the state should be zero
    assert np.allclose(pde.state.vector().get_local(), (end - start - dt) * value)

    # Solving for longer time should not change the state
    pde.solve((end + dt, 2 * end), dt=dt)
    assert np.allclose(pde.state.vector().get_local(), (end - start - dt) * value)


def test_double_stimulation():
    mesh = dolfin.UnitIntervalMesh(10)
    dt = 0.01
    value1 = 2.0
    value2 = 3.0
    start1 = 0.5
    end1 = 1.0
    start2 = 0.9
    end2 = 1.5

    time = dolfin.Constant(0.0)
    expr1 = dolfin.Expression(
        "(t >= start1 && t <= end1) ? value1 : 0.0",
        value1=value1,
        start1=start1,
        end1=end1,
        t=time,
        degree=0,
    )

    expr2 = dolfin.Expression(
        "(t >= start2 && t <= end2) ? value2 : 0.0",
        value2=value2,
        start2=start2,
        end2=end2,
        t=time,
        degree=0,
    )

    I_s = [
        beat.stimulation.Stimulus(dz=dolfin.dx, expr=expr1),
        beat.stimulation.Stimulus(dz=dolfin.dx, expr=expr2),
    ]

    pde = beat.MonodomainModel(
        time=time,
        mesh=mesh,
        M=dolfin.Constant(0.0),
        I_s=I_s,
        params=dict(theta=1.0),
    )

    pde.step((0.0, 0.4))
    assert np.allclose(pde.state.vector().get_local(), 0.0)

    # Solve up to the second stimulus starts
    t0 = 0.9
    stim_t0 = value1 * (t0 - start1)
    pde.solve((0.4, t0), dt=dt)
    # At time dt the stimulus should be value and since M is zero the state should be value * dt
    assert np.allclose(pde.state.vector().get_local(), stim_t0)

    # Solve up to the end of the first stimulus
    pde.solve((t0, end1 + dt), dt=dt)
    assert np.allclose(
        pde.state.vector().get_local(),
        (end1 - start1 - dt) * value1 + (end1 + dt - start2) * value2,
    )

    # Solve up to the end of the second stimulus
    pde.solve((end1 + dt, end2 + dt), dt=dt)
    assert np.allclose(
        pde.state.vector().get_local(),
        (end1 - start1 - dt) * value1 + (end2 - start2 - 2 * dt) * value2,
    )

    # Solving for longer time should not change the state
    pde.solve((end2 + dt, 2 * end2), dt=dt)
    assert np.allclose(
        pde.state.vector().get_local(),
        (end1 - start1 - dt) * value1 + (end2 - start2 - 2 * dt) * value2,
    )
