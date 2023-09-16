from typing import NamedTuple, Callable
from scipy.integrate import solve_ivp
import numpy as np
import numpy.typing as npt

from .cellmodel import CellModel


class ODEResults(NamedTuple):
    y: npt.NDArray[np.float64]
    t: npt.NDArray[np.float64]


def solve_cellmodel(
    model: CellModel,
    interval: tuple[float, float],
    i_app: Callable[[float], float] | None = None,
    **kwargs,
) -> ODEResults:
    if i_app is None:
        i_app = lambda t: 0.0

    def fun(t, y):
        v = -model.I(y[0:1], y[1:], t) + i_app(t)
        s = model.F(y[0:1], y[1:], t)

        return [*v, *s]

    res = solve_ivp(
        fun,
        interval,
        list(model.ic.values()),
        **kwargs,
    )
    return ODEResults(y=res.y, t=res.t)
