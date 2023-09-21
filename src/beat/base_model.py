from __future__ import annotations
from typing import Any, NamedTuple
import abc
import logging
from enum import Enum, auto

import dolfin
import ufl_legacy as ufl
from ufl_legacy.core.expr import Expr


logger = logging.getLogger(__name__)


class Status(str, Enum):
    OK = auto()
    NOT_CONVERGING = auto()


class Results(NamedTuple):
    state: dolfin.Function
    status: Status


class BaseModel:
    def __init__(
        self,
        time: dolfin.Constant,
        mesh: dolfin.Mesh,
        params: dict[str, Any] | None = None,
    ) -> None:
        self._mesh = mesh
        self.time = time

        self.parameters = type(self).default_parameters()
        if params is not None:
            self.parameters.update(params)

        self._setup_state_space()

        self._timestep = dolfin.Constant(self.parameters["default_timestep"])
        (G, self._prec) = self.variational_forms(self._timestep)
        self._lhs, self._rhs = dolfin.system(G)

        logger.debug("Preassembling monodomain matrix (and initializing vector)")
        self._lhs_matrix = dolfin.assemble(self._lhs)
        self._rhs_vector = dolfin.Vector(
            self._mesh.mpi_comm(), self._lhs_matrix.size(0)
        )
        self._lhs_matrix.init_vector(self._rhs_vector, 0)

        # Create linear solver (based on parameter choices)
        self.linear_solver, self._update_solver = self._create_linear_solver()

    @abc.abstractmethod
    def _setup_state_space(self) -> None:
        ...

    def _create_linear_solver(self):
        "Helper function for creating linear solver based on parameters."
        solver_type = self.parameters["linear_solver_type"]

        if solver_type == "direct":
            solver = dolfin.LUSolver(self._lhs_matrix, self.parameters["lu_type"])
            solver.parameters.update(self.parameters["lu_solver"])
            update_routine = self._update_lu_solver

        elif solver_type == "iterative":
            # Preassemble preconditioner (will be updated if time-step
            # changes)
            logger.debug("Preassembling preconditioner")
            # Initialize KrylovSolver with matrix and preconditioner
            alg = self.parameters["algorithm"]
            prec = self.parameters["preconditioner"]
            if self.parameters["use_custom_preconditioner"]:
                self._prec_matrix = dolfin.assemble(self._prec)
                solver = dolfin.PETScKrylovSolver(alg, prec)
                solver.parameters.update(self.parameters["krylov_solver"])
                solver.set_operators(self._lhs_matrix, self._prec_matrix)
                solver.ksp().setFromOptions()
            else:
                solver = dolfin.PETScKrylovSolver(alg, prec)
                solver.parameters.update(self.parameters["krylov_solver"])
                solver.set_operator(self._lhs_matrix)
                solver.ksp().setFromOptions()

            update_routine = self._update_krylov_solver
        else:
            msg = f"Unknown linear_solver_type given: {solver_type}"
            raise ValueError(msg)

        return (solver, update_routine)

    @property
    @abc.abstractmethod
    def state(self) -> dolfin.Function:
        ...

    @abc.abstractmethod
    def assign_previous(self) -> None:
        ...

    @staticmethod
    def default_parameters():
        def to_dict(d):
            if isinstance(d, dolfin.Parameters):
                return to_dict(d.to_dict())

            elif isinstance(d, dict):
                res = {}
                for k, v in d.items():
                    res[k] = to_dict(v)
                return res
            elif isinstance(d, dolfin.cpp.parameter.Parameter):
                return d.value()
            else:
                return d

        lu_solver_params = to_dict(dolfin.LUSolver.default_parameters())
        krylov_solver_params = to_dict(dolfin.KrylovSolver.default_parameters())

        return {
            "theta": 0.5,
            "degree": 1,
            "family": "Lagrange",
            "default_timestep": 1.0,
            "linear_solver_type": "iterative",
            "lu_type": "default",
            "algorithm": "cg",
            "preconditioner": "petsc_amg",
            "krylov_solver": krylov_solver_params,
            "lu_solver": lu_solver_params,
        }

    @abc.abstractmethod
    def variational_forms(self, k_n: Expr | float) -> tuple[ufl.Form, ufl.Form]:
        """Create the variational forms corresponding to the given
        discretization of the given system of equations.

        *Arguments*
          k_n (:py:class:`ufl.Expr` or float)
            The time step

        *Returns*
          (G, prec) (:py:class:`tuple` of :py:class:`ufl.Form`)

        """
        ...

    def step(self, interval):
        """
        Solve on the given time step (t0, t1).

        *Arguments*
          interval (:py:class:`tuple`)
            The time interval (t0, t1) for the step

        *Invariants*
          Assuming that v_ is in the correct state for t0, gives
          self.v in correct state at t1.
        """

        # timer = dolfin.Timer("PDE Step")

        # Extract interval and thus time-step
        (t0, t1) = interval
        dt = t1 - t0
        theta = self.parameters["theta"]
        t = t0 + theta * dt
        self.time.assign(t)

        # Update matrix and linear solvers etc as needed
        timestep_unchanged = abs(dt - float(self._timestep)) < 1.0e-12
        self._update_solver(timestep_unchanged, dt)

        # Assemble right-hand-side
        dolfin.assemble(self._rhs, tensor=self._rhs_vector)

        # Solve problem
        self.linear_solver.solve(self.state.vector(), self._rhs_vector)
        # timer.stop()

    def _update_lu_solver(self, timestep_unchanged, dt):
        """Helper function for updating an LUSolver depending on
        whether timestep has changed."""

        if timestep_unchanged:
            logger.debug("Timestep is unchanged, reusing LU factorization")
        else:
            logger.debug("Timestep has changed, updating LU factorization")
            # Update stored timestep
            # FIXME: dolfin_adjoint still can't annotate constant assignment.
            self._timestep.assign(dolfin.Constant(dt))  # , annotate=annotate)

            # Reassemble matrix
            dolfin.assemble(self._lhs, tensor=self._lhs_matrix)

    def _update_krylov_solver(self, timestep_unchanged, dt):
        """Helper function for updating a KrylovSolver depending on
        whether timestep has changed."""

        # Update reuse of preconditioner parameter in accordance with
        # changes in timestep
        if timestep_unchanged:
            logger.debug("Timestep is unchanged, reusing preconditioner")
            # self.linear_solver.parameters["preconditioner"]["structure"] = "same"
        else:
            logger.debug("Timestep has changed, updating preconditioner")
            # self.linear_solver.parameters["preconditioner"]["structure"] = \
            #                                            "same_nonzero_pattern"

            # Update stored timestep
            self._timestep.assign(dolfin.Constant(dt))

            # Reassemble matrix
            dolfin.assemble(self._lhs, tensor=self._lhs_matrix)

            # Reassemble preconditioner
            if self.parameters["use_custom_preconditioner"]:
                dolfin.assemble(self._prec, tensor=self._prec_matrix)

        if self.state.vector().norm("l2") > 1.0e-12:
            logger.debug("Initial guess is non-zero.")
            self.linear_solver.parameters["nonzero_initial_guess"] = True

    def solve(
        self,
        interval: tuple[float, float],
        dt: float | None = None,
    ):
        """
        Solve the discretization on a given time interval (t0, t1)
        with a given timestep dt and return generator for a tuple of
        the interval and the current solution.

        *Arguments*
          interval (:py:class:`tuple`)
            The time interval for the solve given by (t0, t1)
          dt (int, optional)
            The timestep for the solve. Defaults to length of interval

        *Returns*
          (timestep, solution_field) via (:py:class:`genexpr`)

        *Example of usage*::

          # Create generator
          solutions = solver.solve((0.0, 1.0), 0.1)

          # Iterate over generator (computes solutions as you go)
          for (interval, solution_fields) in solutions:
            (t0, t1) = interval
            v_, v = solution_fields
            # do something with the solutions
        """

        # Initial set-up
        # Solve on entire interval if no interval is given.
        (T0, T) = interval
        if dt is None:
            dt = T - T0
        t0 = T0
        t1 = T0 + dt

        # Step through time steps until at end time
        while True:
            logger.info("Solving on t = (%g, %g)" % (t0, t1))
            self.step((t0, t1))

            # Yield solutions
            # yield (t0, t1), self.solution_fields()

            # Break if this is the last step
            if (t1 + dt) > (T + dolfin.DOLFIN_EPS):
                break

            self.assign_previous()

            t0 = t1
            t1 = t0 + dt

        return Results(state=self.state, status=Status.OK)
