from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple, Callable

import numpy as np
import numpy.typing as npt
import dolfin

from .utils import local_project

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
class ODESystemSolver:
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
        self.states[:] = self.fun(
            states=self.states, t=t0, parameters=self.parameters, dt=dt
        )


@dataclass
class DolfinODESolver:
    v_ode: dolfin.Function
    v_pde: dolfin.Function
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

        self._ode = ODESystemSolver(
            fun=self.fun,
            states=self._values,
            parameters=self.parameters,
        )

    def to_dolfin(self) -> None:
        """Assign values from numpy array to dolfin function"""
        self.v_ode.vector().set_local(self._values[self.v_index, :])

    def from_dolfin(self) -> None:
        """Assign values from dolfin function to numpy array"""
        self._values[self.v_index, :] = self.v_ode.vector().get_local()

    # ode_to_pde projects v_ode (DG0, quadrature space, ...) into v_pde (CG1)
    def ode_to_pde(self) -> None:
        local_project(self.v_ode, self.v_pde.function_space(), self.v_pde)

    # pde_to_ode projects v_pde (CG1) into v_ode (DG0, quadrature space, ...)
    def pde_to_ode(self) -> None:
        local_project(self.v_pde, self.v_ode.function_space(), self.v_ode)

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
        return self.v_ode.vector().size()

    def step(self, t0: float, dt: float):
        self._ode.step(t0=t0, dt=dt)

    @property
    def full_values(self):
        return self._values


@dataclass
class DolfinMultiODESolver:
    v_ode: dolfin.Function
    v_pde: dolfin.Function
    markers: dolfin.Function
    init_states: dict[int, npt.NDArray]
    parameters: dict[int, npt.NDArray]
    fun: dict[int, Callable]
    num_states: dict[int, int]
    v_index: dict[int, int]

    def __post_init__(self):
        if self.v_ode.vector().size() != self.markers.vector().size():
            raise RuntimeError(
                "Marker and voltage need to be in the same function space"
            )

        self._marker_values = tuple(self.init_states.keys())
        self._num_points = {}
        self._odes = {}
        self._values = {}
        self._inds = {}

        self._initialize_full_values()

        for marker in self._marker_values:
            where = self.markers.vector().get_local() == marker
            self._num_points[marker] = where.sum()
            self._inds[marker] = where

            if np.shape(self.init_states[marker]) == self.shape(marker):
                self._values[marker] = np.copy(self.init_states[marker])
            else:
                self._values[marker] = np.zeros(self.shape(marker))
                self._values[marker].T[:] = self.init_states[marker]

            self._odes[marker] = ODESystemSolver(
                fun=self.fun[marker],
                states=self._values[marker],
                parameters=self.parameters[marker],
            )

    def _initialize_full_values(self):
        self._all_states_equal_size = (
            np.array(tuple(self.num_states.values()))
            == tuple(self.num_states.values())[0]
        ).all()
        if self._all_states_equal_size:
            self._full_values = np.zeros(
                (next(iter(self.num_states.values())), self.markers.vector().size())
            )

    def to_dolfin(self) -> None:
        """Assign values from numpy array to dolfin function"""
        arr = self.v_ode.vector().get_local().copy()
        for marker in self._marker_values:
            arr[self._inds[marker]] = self._values[marker][self.v_index[marker], :]
        self.v_ode.vector().set_local(arr)

    def from_dolfin(self) -> None:
        """Assign values from dolifn function to numpy array"""
        arr = self.v_ode.vector().get_local()
        for marker in self._marker_values:
            self._values[marker][self.v_index[marker], :] = arr[self._inds[marker]]

    # ode_to_pde projects v_ode (DG0, quadrature space, ...) into v_pde (CG1)
    def ode_to_pde(self) -> None:
        local_project(self.v_ode, self.v_pde.function_space(), self.v_pde)

    # pde_to_ode projects v_pde (CG1) into v_ode (DG0, quadrature space, ...)
    def pde_to_ode(self) -> None:
        local_project(self.v_pde, self.v_ode.function_space(), self.v_ode)

    def values(self, marker: int) -> npt.NDArray:
        return self._values[marker]

    def num_parameters(self, marker: int) -> int:
        return len(self.parameters[marker])

    def shape(self, marker: int) -> tuple[int, int]:
        return (self.num_states[marker], self._num_points[marker])

    def num_points(self, marker: int) -> int:
        return self._num_points[marker]

    def step(self, t0: float, dt: float):
        for ode in self._odes.values():
            ode.step(t0=t0, dt=dt)

    @property
    def full_values(self):
        if not self._all_states_equal_size:
            msg = (
                "Cannot get full values size states are not of equal size. "
                f"Have {self.num_states=}, use .values(marker) instead"
            )
            raise RuntimeError(msg)

        for marker in self._marker_values:
            where = self.markers.vector().get_local() == marker
            self._full_values[:, where] = self._values[marker]

        return self._full_values
