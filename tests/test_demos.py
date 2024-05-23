"""Run all demos and mock all solve calls"""

import pytest
import contextlib
from pathlib import Path
import os
from unittest.mock import patch
import sys


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


@pytest.mark.parametrize(
    "demo",
    [
        "diffusion.py",
        "conduction_velocty.py",
        "endocardial_stimulation.py",
        "niederer_benchmark.py",
    ],
)
def test_monodomain_demos(demo):
    with chdir(DEMODIR):
        print(os.getcwd())
        with patch("beat.monodomain_solver.MonodomainSplittingSolver.step"):
            __import__(demo[:-3])
