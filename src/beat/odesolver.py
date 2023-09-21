from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple, Callable

import numpy as np
import numpy.typing as npt
import dolfin


EPS = 1e-12


class ODEResults(NamedTuple):
    y: npt.NDArray[np.float64]
    t: npt.NDArray[np.float64]


def solve(
    fun: np.NDArray,
    t_bound: float,
    states: np.NDArray,
    V: np.NDArray,
    V_index: int,
    dt: float,
    parameters: np.NDArray,
    t0: float = 0.0,
    extra: dict[str, float | npt.NDArray] | None = None,
):
    if extra is None:
        extra = {}
    i = 0
    t = t0
    while t + dt < t_bound:
        fun(states=states, t=t, parameters=parameters, dt=dt, **extra)
        V[i, :] = states[V_index, :]
        i += 1
        t += dt


@dataclass
class ODESytemSolver:
    fun: Callable
    states: npt.NDArray
    dt: float
    t_bound: float
    parameters: npt.NDArray
    t0: float = 0

    def __post_init__(self):
        self.t = self.t0

    @property
    def num_points(self) -> int:
        return self.states.shape[1]

    @property
    def num_states(self) -> int:
        return self.states.shape[0]

    def step(self) -> None:
        if self.t + self.dt >= self.t_bound + EPS:
            return
        self.fun(states=self.states, t=self.t, parameters=self.parameters, dt=self.dt)
        self.t += self.dt


class Assigner(NamedTuple):
    f: dolfin.Function
    assigner_to_single: dolfin.FunctionAssigner
    assigner_from_single: dolfin.FunctionAssigner


def setup_assigner(vs, index):
    # Set-up separate potential function for post processing
    VS0 = vs.function_space().sub(index)
    V = VS0.collapse()
    v = dolfin.Function(V)
    # Set-up object to optimize assignment from a function to subfunction
    assigner_to_single = dolfin.FunctionAssigner(V, VS0)
    assigner_to_single.assign(v, vs.sub(index))
    assigner_from_single = dolfin.FunctionAssigner(VS0, V)
    assigner_from_single.assign(vs.sub(index), v)
    return Assigner(
        f=v,
        assigner_to_single=assigner_to_single,
        assigner_from_single=assigner_from_single,
    )


@dataclass
class DolfinODESolver:
    s: dolfin.Function
    init_states: npt.NDArray
    parameters: npt.NDArray
    fun: Callable
    dt: float
    t_bound: float
    t0: float = 0.0

    def __post_init__(self):
        self.assigners = [
            setup_assigner(self.s, index) for index in range(self.num_states)
        ]
        self._values = np.zeros((self.num_states, self.num_points))
        self._values.T[:] = self.init_states
        self._ode = ODESytemSolver(
            fun=self.fun,
            states=self._values,
            dt=self.dt,
            t_bound=self.t_bound,
            t0=self.t0,
            parameters=self.parameters,
        )

    def to_dolfin(self):
        """Assign values from numpy array to dolfin function"""
        for i, a in enumerate(self.assigners):
            a.f.vector().set_local(self._values[i, :])
            a.assigner_from_single.assign(self.s.sub(i), a.f)

    @property
    def values(self):
        return self._values

    @property
    def num_states(self) -> int:
        return self.s.num_sub_spaces()

    @property
    def t(self) -> float:
        return self._ode.t

    @property
    def num_points(self) -> int:
        # FIXME: Should be number of dofs in order to work with MPI
        return self.s.function_space().sub(0).collapse().dim()

    def step(self):
        self._ode.step()
