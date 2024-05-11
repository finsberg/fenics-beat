from pathlib import Path
from beat.single_cell import get_steady_state, compute_hash
import gotranx
import numpy as np

import pytest

here = Path(__file__).resolve().parent


@pytest.fixture
def model():
    odefile = here / ".." / "odes" / "ORdmm_Land.ode"
    ode = gotranx.load_ode(odefile)
    code = gotranx.cli.gotran2py.get_code(
        ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
    )
    model = {}
    exec(code, model)
    return model


def test_get_steady_ORd_with_tracked_values(tmp_path, model):
    track_indices = [model["state_index"]("v"), model["state_index"]("cai")]

    nbeats = 2
    BCL = 1000
    dt = 0.05
    save_every_ms = 2.0

    outdir = tmp_path / "endo"
    y = get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=0),
        outdir=outdir,
        track_indices=track_indices,
        nbeats=nbeats,
        save_every_ms=save_every_ms,
        BCL=BCL,
        dt=dt,
    )
    hash_input = compute_hash(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=0),
        nbeats=nbeats,
        BCL=BCL,
        dt=dt,
    )

    assert (outdir / f"steady_states_{hash_input}.npy").is_file()
    assert (outdir / f"tracked_values_{hash_input}.npy").is_file()
    assert (outdir / f"tracked_values_{hash_input}.png").is_file()

    tracked_values = np.load(outdir / f"tracked_values_{hash_input}.npy")
    assert tracked_values.shape == (
        round(nbeats * BCL // save_every_ms),
        len(track_indices),
    )
    assert np.isclose(tracked_values[0, 0], model["init_state_values"]()[track_indices[0]])
    assert not np.isclose(tracked_values[0, 0], tracked_values[1, 0])

    y = np.load(outdir / f"steady_states_{hash_input}.npy")
    assert y.shape == (len(model["init_state_values"]()),)
    assert not np.allclose(y, model["init_state_values"]())


def test_get_steady_ORd_without_tracked_values(tmp_path, model):
    nbeats = 2
    BCL = 1000
    dt = 0.05

    outdir = tmp_path / "epi"
    y = get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=0),
        outdir=outdir,
        nbeats=nbeats,
        BCL=BCL,
        dt=dt,
    )
    hash_input = compute_hash(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=model["init_parameter_values"](celltype=0),
        nbeats=nbeats,
        BCL=BCL,
        dt=dt,
    )

    assert (outdir / f"steady_states_{hash_input}.npy").is_file()
    assert not (outdir / f"tracked_values_{hash_input}.npy").is_file()
    assert not (outdir / f"tracked_values_{hash_input}.png").is_file()

    y = np.load(outdir / f"steady_states_{hash_input}.npy")
    assert y.shape == (len(model["init_state_values"]()),)
    assert not np.allclose(y, model["init_state_values"]())
