import logging
import ufl_legacy as ufl
import dolfin

from .base_model import BaseModel


logger = logging.getLogger(__name__)


class BidomainModel(BaseModel):
    def __init__(self, mesh, time, M_i, M_e, I_s=None, I_a=None, v_=None, params=None):
        super().__init__(mesh=mesh, params=params)

        self._nullspace_basis = None
        # Store input

        self.time = time
        self._M_i = M_i
        self._M_e = M_e
        self._I_s = I_s
        self._I_a = I_a

        # Set-up function spaces
        k = self.parameters["degree"]
        Ve = dolfin.FiniteElement("CG", self._mesh.ufl_cell(), k)
        V = dolfin.FunctionSpace(self._mesh, "CG", k)
        Ue = dolfin.FiniteElement("CG", self._mesh.ufl_cell(), k)
        dolfin.FunctionSpace(self._mesh, "CG", k)

        use_R = self.parameters["use_avg_u_constraint"]
        if use_R:
            Re = dolfin.FiniteElement("R", self._mesh.ufl_cell(), 0)
            dolfin.FunctionSpace(self._mesh, "R", 0)
            self.VUR = dolfin.FunctionSpace(mesh, dolfin.MixedElement((Ve, Ue, Re)))
        else:
            self.VUR = dolfin.FunctionSpace(mesh, dolfin.MixedElement((Ve, Ue)))

        self.V = V

        # Set-up solution fields:
        self.merger = dolfin.FunctionAssigner(V, self.VUR.sub(0))
        self.v_ = dolfin.Function(V, name="v_")

        self._state = dolfin.Function(self.VUR, name="vur")

    @property
    def state(self) -> dolfin.Function:
        return self._state

    def assign_previous(self) -> None:
        self.merger.assign(self.v_, self.state.sub(0))

    def _create_linear_solver(self):
        "Helper function for creating linear solver based on parameters."
        solver, update_routine = super()._create_linear_solver()
        solver_type = self.parameters["linear_solver_type"]

        if solver_type == "iterative":
            # Still waiting for that bug fix:
            solver.parameters.update({"convergence_norm_type": "preconditioned"})
            solver.parameters.update(self.parameters["petsc_krylov_solver"])

            # Set nullspace if present. We happen to know that the
            # transpose nullspace is the same as the nullspace (easy
            # to prove from matrix structure).
            if not self.parameters["use_avg_u_constraint"]:
                A = dolfin.as_dolfin_type(self._lhs_matrix)
                A.set_nullspace(self.nullspace)
        return solver, update_routine

    @property
    def nullspace(self):
        if self._nullspace_basis is None:
            null_vector = dolfin.Vector(self.state.vector())
            self.VUR.sub(1).dofmap().set(null_vector, 1.0)
            null_vector *= 1.0 / null_vector.norm("l2")
            self._nullspace_basis = dolfin.VectorSpaceBasis([null_vector])
        return self._nullspace_basis

    @staticmethod
    def default_parameters():
        """Initialize and return a set of default parameters

        *Returns*
          A set of parameters (:py:class:`dolfin.Parameters`)

        To inspect all the default parameters, do::

          info(BidomainSolver.default_parameters(), True)
        """
        params = super(BidomainModel, BidomainModel).default_parameters()
        params["use_avg_u_constraint"] = False
        return params

    def variational_forms(
        self, k_n: ufl.core.expr.Expr | float
    ) -> tuple[ufl.Form, ufl.Form]:
        """Create the variational forms corresponding to the given
        discretization of the given system of equations.

        *Arguments*
          k_n (:py:class:`ufl.Expr` or float)
            The time step

        *Returns*
          (lhs, rhs) (:py:class:`tuple` of :py:class:`ufl.Form`)

        """

        # Extract theta parameter and conductivities
        theta = self.parameters["theta"]
        M_i = self._M_i
        M_e = self._M_e

        # Define variational formulation
        use_R = self.parameters["use_avg_u_constraint"]
        if use_R:
            (v, u, l) = dolfin.TrialFunctions(self.VUR)
            (w, q, lamda) = dolfin.TestFunctions(self.VUR)
        else:
            (v, u) = dolfin.TrialFunctions(self.VUR)
            (w, q) = dolfin.TestFunctions(self.VUR)

        # Set-up measure and rhs from stimulus
        # (dz, rhs) = rhs_with_markerwise_field(self._I_s, self._mesh, w)
        dz = dolfin.dx
        rhs = self._I_s * w * dz()

        # Set-up variational problem
        Dt_v_k_n = v - self.v_
        v_mid = theta * v + (1.0 - theta) * self.v_
        theta_parabolic = (
            ufl.inner(M_i * ufl.grad(v_mid), ufl.grad(w)) * dz()
            + ufl.inner(M_i * ufl.grad(u), ufl.grad(w)) * dz()
        )
        theta_elliptic = (
            ufl.inner(M_i * ufl.grad(v_mid), ufl.grad(q)) * dz()
            + ufl.inner((M_i + M_e) * ufl.grad(u), ufl.grad(q)) * dz()
        )
        G = (
            Dt_v_k_n * w * dz()
            + k_n * theta_parabolic
            + k_n * theta_elliptic
            - k_n * rhs
        )

        if use_R:
            G += k_n * (lamda * u + l * q) * dz()

        # Add applied current as source in elliptic equation if
        # applicable
        if self._I_a:
            G -= k_n * self._I_a * q * dz()

        # (a, L) = dolfin.system(G)
        # return (a, L)
        return G, None
