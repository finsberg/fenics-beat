import numpy as np

import beat


def test_fitzhughnagumo():
    model = beat.cellmodels.FitzHughNagumo()

    i_app = lambda t: 0.05 * 125 if 50 <= t <= 60 else 0
    x = np.linspace(0, 1000, 100)
    res = beat.cellsolver.solve_cellmodel(
        model, [0, 1000], i_app=i_app, t_eval=x, max_step=1
    )

    assert np.isclose(res.y.min(), -85)
    assert np.isclose(res.y.max(), 68.60182926019024)
    assert np.isclose(res.t[0], 0.0)
    assert np.isclose(res.t[-1], 1000.0)
    # import matplotlib.pyplot as plt
    # plt.plot(res.t, res.y[0, :])
    # plt.savefig("cell.png")
