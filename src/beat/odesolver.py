from __future__ import annotations
from dataclasses import dataclass
from typing import NamedTuple, Callable

import numpy as np
import numpy.typing as npt

EPS = 1e-12


class ODEResults(NamedTuple):
    y: npt.NDArray[np.float64]
    t: npt.NDArray[np.float64]


# def solve_cellmodel(
#     model: CellModel,
#     interval: tuple[float, float],
#     y0: npt.NDArray[float],
#     i_app: Callable[[float], float] | None = None,
#     method="Radau",
#     **kwargs,
# ) -> ODEResults:
#     if i_app is None:
#         i_app = lambda t: 0.0

#     def fun(t, y):
#         v = -model.I(y[:, 0:1], y[:, 1:], t) + i_app(t)
#         s = model.F(y[:, 0:1], y[:, 1:], t)
#         return np.hstack((v, s))

#     # res = solve_ivp(
#     #     fun,
#     #     interval,
#     #     list(model.y0.values()),
#     #     **kwargs,
#     # )
#     t0, t_bound = interval
#     if method == "RK45":
#         ode = RK45(fun, t0=t0, y0=y0, t_bound=t_bound, **kwargs)
#     else:
#         ode = Radau(fun, t0=t0, y0=y0, t_bound=t_bound, **kwargs)

#     status = None
#     ts, ys = [], []
#     while status is None:
#         message = ode.step()
#         if ode.status == "finished":
#             status = 0
#         elif ode.status == "failed":
#             status = -1
#             break

#         ts.append(ode.t)
#         ys.append(ode.y)

#     return ODEResults(y=np.array(ys), t=np.array(ts))


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
