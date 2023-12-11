import logging
import dolfin
import numpy.typing as npt
import numpy as np
import ufl_legacy as ufl


logger = logging.getLogger(__name__)


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
