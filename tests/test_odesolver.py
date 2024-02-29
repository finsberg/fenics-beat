import numpy as np
import dolfin

import beat
from beat.odesolver import ODESystemSolver, DolfinODESolver, DolfinMultiODESolver


def test_simple_ode_odesystemsolver():
    def simple_ode_forward_euler(states, t, dt, parameters):
        v, s = states
        values = np.zeros_like(states)
        values[0] = v - s * dt
        values[1] = s + v * dt
        return values

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
        ode = ODESystemSolver(
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
    model = beat.cellmodels.beeler_reuter_1977
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

    ode = ODESystemSolver(
        fun=beat.cellmodels.beeler_reuter_1977.forward_generalized_rush_larsen,
        states=states,
        parameters=parameters,
    )
    assert np.allclose(ode.states, old_states)

    ode.step(t0, dt)

    assert not np.allclose(ode.states, old_states)


def test_beeler_reuter_unit_square():
    model = beat.cellmodels.beeler_reuter_1977
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
    parameters[model.parameter_indices("IstimAmplitude")] = 1.0

    mesh = dolfin.UnitSquareMesh(5, 5)
    V = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    s = dolfin.Function(V)
    dt = 0.1
    t0 = 0.0

    dolfin_ode = DolfinODESolver(
        v_ode=dolfin.Function(V),
        v_pde=s,
        num_states=len(init_states),
        fun=beat.cellmodels.beeler_reuter_1977.forward_generalized_rush_larsen,
        init_states=init_states,
        parameters=parameters,
    )
    assert np.allclose(dolfin_ode.v_ode.vector().get_local(), 0.0)
    dolfin_ode.to_dolfin()
    dolfin_ode.ode_to_pde()

    # Just check that values have been updated
    old_state = dolfin_ode.v_pde.vector().get_local().copy()
    assert not np.allclose(old_state, 0.0)

    N = 10
    t = t0
    for _ in range(N):
        dolfin_ode.step(t, dt)
        t += dt

    dolfin_ode.to_dolfin()
    dolfin_ode.ode_to_pde()
    assert not np.allclose(dolfin_ode.v_pde.vector().get_local(), old_state)


def test_assignment_ode():
    model = beat.cellmodels.beeler_reuter_1977
    init_states = model.init_state_values()
    parameters = model.init_parameter_values()
    parameters[model.parameter_indices("IstimAmplitude")] = 1.0
    v_index = model.state_indices("V")

    mesh = dolfin.UnitSquareMesh(5, 5)
    V = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    v = dolfin.Function(V)
    ode = DolfinODESolver(
        v_ode=dolfin.Function(V),
        v_pde=v,
        num_states=len(init_states),
        fun=beat.cellmodels.beeler_reuter_1977.forward_generalized_rush_larsen,
        init_states=init_states,
        parameters=parameters,
        v_index=v_index,
    )
    assert np.allclose(ode.v_ode.vector().get_local(), 0)
    assert np.allclose(ode.values[:, 0], init_states)

    ode.to_dolfin()
    ode.ode_to_pde()
    assert np.allclose(ode.v_pde.vector().get_local(), init_states[v_index])

    # Now update values for v
    ode.values[v_index, :] = 42.0
    assert np.allclose(ode.v_ode.vector().get_local(), init_states[v_index])
    ode.to_dolfin()
    ode.ode_to_pde()
    assert np.allclose(ode.v_pde.vector().get_local(), 42.0)

    # Now update dolfin function for v
    ode.v_pde.assign(dolfin.Constant(13.0))

    ode.pde_to_ode()
    ode.from_dolfin()

    assert np.allclose(ode.values[v_index, :], 13.0)
    assert np.allclose(ode.full_values[v_index, :], 13.0)


def test_ode_with_markers_3D_to_and_from_dolfin():
    model = beat.cellmodels.tentusscher_panfilov_2006
    mesh = dolfin.UnitCubeMesh(3, 3, 3)
    V = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    markers = dolfin.Function(V)
    arr = markers.vector().get_local().copy()
    v2d = dolfin.vertex_to_dof_map(V)

    mfun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    mfun.set_all(0)
    dolfin.CompiledSubDomain("near(x[0], 0)").mark(mfun, 1)
    dolfin.CompiledSubDomain("near(x[0], 1)").mark(mfun, 2)

    for marker in [0, 1, 2]:
        for facet in mfun.where_equal(marker):
            f = dolfin.Facet(mesh, facet)
            arr[v2d[f.entities(0)]] = marker
    markers.vector().set_local(arr)

    v = dolfin.Function(V)
    init_states = {
        0: model.mid.init_state_values(),
        1: model.endo.init_state_values(),
        2: model.epi.init_state_values(),
    }
    parameters = {
        0: model.mid.init_parameter_values(),
        1: model.endo.init_parameter_values(),
        2: model.epi.init_parameter_values(),
    }
    fun = {
        0: model.mid.forward_generalized_rush_larsen,
        1: model.endo.forward_generalized_rush_larsen,
        2: model.epi.forward_generalized_rush_larsen,
    }
    v_index = {
        0: model.mid.state_indices("V"),
        1: model.endo.state_indices("V"),
        2: model.epi.state_indices("V"),
    }

    ode = DolfinMultiODESolver(
        v_ode=dolfin.Function(V),
        v_pde=v,
        markers=markers,
        num_states={i: len(s) for i, s in init_states.items()},
        fun=fun,
        init_states=init_states,
        parameters=parameters,
        v_index=v_index,
    )

    marker_values = {0: 40, 1: 41, 2: 42}

    for m, v in marker_values.items():
        ode.values(m)[ode.v_index[m], :] = v
    ode.to_dolfin()
    ode.ode_to_pde()

    for m, v in marker_values.items():
        assert np.allclose(ode.v_pde.vector()[markers.vector().get_local() == m], v)

    # Now go the other way
    marker_values = {0: 4, 1: 5, 2: 6}
    arr = ode.v_pde.vector().get_local().copy()
    for m, v in marker_values.items():
        arr[markers.vector().get_local() == m] = v
    ode.v_pde.vector().set_local(arr)
    ode.pde_to_ode()
    ode.from_dolfin()

    v_arr = np.zeros_like(arr)
    for m, v in marker_values.items():
        assert np.allclose(ode.values(m)[v_index[m], :], v)
        v_arr[markers.vector().get_local() == m] = v

    assert ode.full_values.shape == (len(init_states[0]), len(arr))
    assert np.allclose(ode.full_values[v_index[0], :], v_arr)
