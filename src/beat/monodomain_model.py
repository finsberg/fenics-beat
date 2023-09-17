import logging

import dolfin
import ufl_legacy as ufl
from ufl_legacy.core.expr import Expr

from .base_model import BaseModel

logger = logging.getLogger(__name__)


class MonodomainModel(BaseModel):
    r"""Solve

    .. math::

        \chi C_m \frac{\partial V}{\partial t} -
        \nabla \cdot \left( M \nabla V \right) - I_{\mathrm{stim}} = 0

    """

    def __init__(
        self, time, mesh: dolfin.Mesh, M: ufl.Coefficient | float, I_s, params=None
    ) -> None:
        self._M = M
        self._I_s = I_s
        self.time = time
        super().__init__(mesh=mesh, params=params)

    def _setup_state_space(self) -> None:
        # Set-up function spaces
        k = self.parameters["degree"]
        family = self.parameters["family"]

        element = dolfin.FiniteElement(
            family=family,
            cell=self._mesh.ufl_cell(),
            degree=k,
            quad_scheme="default",
        )

        self.V = dolfin.FunctionSpace(self._mesh, element)

        # Set-up solution fields:
        self.v_ = dolfin.Function(self.V, name="v_")
        self._state = dolfin.Function(self.V, name="v")

    @property
    def state(self) -> dolfin.Function:
        return self._state

    def assign_previous(self):
        self.v_.assign(self.state)

    @staticmethod
    def default_parameters():
        params = super(MonodomainModel, MonodomainModel).default_parameters()
        params["use_custom_preconditioner"] = True
        return params

    def variational_forms(self, k_n: Expr | float) -> tuple[ufl.Form, ufl.Form]:
        """Create the variational forms corresponding to the given
        discretization of the given system of equations.

        *Arguments*
          k_n (:py:class:`ufl.Expr` or float)
            The time step

        *Returns*
          (lhs, rhs, prec) (:py:class:`tuple` of :py:class:`ufl.Form`)

        """
        theta = self.parameters["theta"]

        # Define variational formulation
        v = dolfin.TrialFunction(self.V)
        w = dolfin.TestFunction(self.V)

        # Set-up variational problem
        Dt_v_k_n = v - self.v_
        v_mid = theta * v + (1.0 - theta) * self.v_

        theta_parabolic = ufl.inner(self._M * ufl.grad(v_mid), ufl.grad(w))

        G = (Dt_v_k_n * w + k_n * theta_parabolic - k_n * self._I_s * w) * dolfin.dx(
            domain=self._mesh
        )

        # Define preconditioner based on educated(?) guess by Marie
        prec = (
            v * w + k_n / 2.0 * ufl.inner(self._M * ufl.grad(v), ufl.grad(w))
        ) * dolfin.dx(domain=self._mesh)

        return (G, prec)
