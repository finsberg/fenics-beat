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
    parameters: npt.NDArray

    @property
    def num_points(self) -> int:
        return self.states.shape[1]

    @property
    def num_states(self) -> int:
        return self.states.shape[0]

    def step(self, t0: float, dt: float) -> None:
        self.fun(states=self.states, t=t0, parameters=self.parameters, dt=dt)


@dataclass
class DolfinODESolver:
    v: dolfin.Function
    init_states: npt.NDArray
    parameters: npt.NDArray
    fun: Callable
    num_states: int
    v_index: int = 0

    def __post_init__(self):
        if np.shape(self.init_states) == self.shape:
            self._values = np.copy(self.init_states)
        else:
            self._values = np.zeros(self.shape)
            self._values.T[:] = self.init_states

        self._ode = ODESytemSolver(
            fun=self.fun,
            states=self._values,
            parameters=self.parameters,
        )

    def to_dolfin(self) -> None:
        """Assign values from numpy array to dolfin function"""
        self.v.vector().set_local(self._values[self.v_index, :])

    def from_dolfin(self) -> None:
        self.values[self.v_index, :] = self.v.vector().get_local()

    @property
    def values(self):
        return self._values

    @property
    def num_parameters(self) -> int:
        return len(self.parameters)

    @property
    def shape(self) -> tuple[int, int]:
        return (self.num_states, self.num_points)

    @property
    def num_points(self) -> int:
        # FIXME: Should be number of dofs in order to work with MPI
        return self.v.vector().size()

    def step(self, t0: float, dt: float):
        self._ode.step(t0=t0, dt=dt)
