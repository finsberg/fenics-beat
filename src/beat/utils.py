from __future__ import annotations
from typing import Any
import logging
from contextlib import contextmanager
import dolfin
import numpy.typing as npt
import numpy as np
from mpi4py import MPI as pyMPI

try:
    import ufl_legacy as ufl
except ImportError:
    import ufl


logger = logging.getLogger(__name__)


def mpi4py_comm(comm):
    """Get mpi4py communicator"""
    try:
        return comm.tompi4py()
    except AttributeError:
        return comm


def peval(f, *x):
    """Parallel synced eval"""
    try:
        yloc = f(*x)
    except RuntimeError:
        yloc = np.inf * np.ones(f.value_shape())

    comm = mpi4py_comm(f.function_space().mesh().mpi_comm())
    yglob = np.zeros_like(yloc)
    comm.Allreduce(yloc, yglob, op=pyMPI.MIN)
    return yglob


def expand_layer(
    markers: dolfin.Function,
    mfun: dolfin.MeshFunction,
    endo_markers: list[int],
    epi_markers: list[int],
    endo_marker: int,
    epi_marker: int,
    endo_size: float,
    epi_size: float,
) -> npt.NDArray:
    """Expand the endo and epi markers to the rest of the mesh
    with a given size

    Parameters
    ----------
    markers : dolfin.Function
        Function where the markers are stored
    mfun : dolfin.MeshFunction
        Mesh function where the markers are stored
    endo_markers: list[int]
        List of markers for the endocardium. Should
        be a subset of the markers in mfun
    epimarkers: list[int]
        List of markers for the epicardium. Should
        be a subset of the markers in mfun
    endo_marker : int
        Marker for the endocardium. Used to set the
        marker on the array that is returned.
    epi_marker : int
        Marker for the epicardium. Used to set the
        marker on the array that is returned.
    endo_size : float
        Size of the endocardium
    epi_size : float
        Size of the epicardium

    Returns
    -------
    npt.NDArray
        Array with the markers
    """
    # Find the rest of the laplace solutions
    V = markers.function_space()
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    arr = markers.vector().get_local().copy()
    sol = dolfin.Function(V)

    a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = v * dolfin.Constant(0) * ufl.dx

    solver = dolfin.PETScKrylovSolver()
    dolfin.PETScOptions.set("ksp_type", "cg")
    dolfin.PETScOptions.set("ksp_norm_type", "unpreconditioned")
    dolfin.PETScOptions.set("ksp_atol", 1e-15)
    dolfin.PETScOptions.set("ksp_rtol", 1e-10)
    dolfin.PETScOptions.set("ksp_max_it", 10_000)
    dolfin.PETScOptions.set("ksp_error_if_not_converged", False)
    # if ksp_monitor:
    #     dolfin.PETScOptions.set("ksp_monitor")
    # if ksp_view:
    #     dolfin.PETScOptions.set("ksp_view")
    dolfin.PETScOptions.set("pc_type", "hypre")
    dolfin.PETScOptions.set("pc_hypre_type", "boomeramg")
    # if pc_view:
    #     dolfin.PETScOptions.set("pc_view")
    solver.set_from_options()
    sol_arrs = []

    for endo in endo_markers:
        for epi in epi_markers:
            sol.vector()[:] = 0

            # Iterate over the three different cases
            logger.info("Solving Laplace equation")
            bcs = [
                dolfin.DirichletBC(V, 0, mfun, endo, "topological"),
                dolfin.DirichletBC(V, 1, mfun, epi, "topological"),
            ]
            A, b = dolfin.assemble_system(a, L, bcs)

            solver.set_operator(A)
            solver.solve(sol.vector(), b)

            sol_arrs.append(sol.vector().get_local())

    # In BiV we have have one epi and two endo solutions
    # We take the minimum of the two endo solutions
    sol_arr = np.min(sol_arrs, axis=0)
    arr[sol_arr < endo_size] = endo_marker
    arr[sol_arr > 1 - epi_size] = epi_marker
    return arr


@contextmanager
def form_compiler_representation(name: str):
    """Context manager to set the form compiler representation"""
    old = dolfin.parameters["form_compiler"]["representation"]
    dolfin.parameters["form_compiler"]["representation"] = name
    yield
    dolfin.parameters["form_compiler"]["representation"] = old


def local_project(
    v: dolfin.Function,
    V: dolfin.FunctionSpace,
    u: dolfin.Function | None = None,
    metadata: dict[str, Any] | None = None,
) -> dolfin.Function | None:
    """Element-wise projection using LocalSolver

    Parameters
    ----------
    v : dolfin.Function
        Function to be projected
    V : dolfin.FunctionSpace
        Function space to project into
    u : dolfin.Function | None, optional
        Optional function to save the projected function, by default None

    Returns
    -------
    dolfin.Function | None
        The projected function
    """

    dv = dolfin.TrialFunction(V)
    v_ = dolfin.TestFunction(V)
    a_proj = ufl.inner(dv, v_) * ufl.dx(metadata=metadata)
    b_proj = ufl.inner(v, v_) * ufl.dx(metadata=metadata)

    if metadata is None:
        solver = dolfin.LocalSolver(a_proj, b_proj)
    else:
        with form_compiler_representation("quadrature"):
            solver = dolfin.LocalSolver(a_proj, b_proj)

    solver.factorize()
    if u is None:
        u = dolfin.Function(V)
    solver.solve_local_rhs(u)
    return u
