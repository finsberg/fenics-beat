"""Run all demos and mock all solve calls"""

import pytest
import contextlib
from pathlib import Path
import os
from unittest.mock import patch
import sys
import shutil


here = Path(__file__).resolve().parent


@contextlib.contextmanager
def chdir(path):
    cwd = os.getcwd()
    os.chdir(str(path))
    sys.path.insert(0, str(path))

    yield
    os.chdir(cwd)
    sys.path.remove(str(path))


DEMODIR = here / ".." / "demos"
ODEDIR = here / ".." / "odes"


@pytest.fixture(scope="session")
def tmpdir_demo(tmpdir_factory):
    tmp = tmpdir_factory.mktemp("src")
    tmp_demo = (tmp / "demos").mkdir()
    shutil.copytree(ODEDIR, tmp / "odes")

    return tmp_demo


def solve_with_save_side_effect(*args, **kwargs):
    y = kwargs["y"]
    track_values = kwargs["track_values"]
    return y, track_values


def solve_without_save_side_effect(*args, **kwargs):
    y = kwargs["y"]
    return y


@pytest.mark.parametrize(
    "demo",
    [
        "diffusion.py",
        "conduction_velocty.py",
        "endocardial_stimulation.py",
        "niederer_benchmark.py",
        "pvc.py",
        "reentry.py",
        "s1s2_tissue.py",
        "apd_restitution.py",
        "multiple_stimulation_sites.py",
        "simple_ode.py",
    ],
)
def test_monodomain_demos(demo, tmpdir_demo):
    # Copy the demo to a temporary directory

    shutil.copy(DEMODIR / demo, tmpdir_demo)

    with chdir(tmpdir_demo):
        with (
            patch("pyvista"),
            patch("beat.monodomain_solver.MonodomainSplittingSolver.step"),
            patch("ap_features.apd") as apd,
            patch("beat.single_cell.solve_with_save") as solve_with_save,
            patch("beat.single_cell.solve_without_save") as solve_without_save,
        ):
            apd.return_value = 100.0
            solve_with_save.side_effect = solve_with_save_side_effect
            solve_without_save.side_effect = solve_without_save_side_effect

            __import__(demo[:-3])
