import numpy as np
import beat


def simple_ode_forward_euler(states, t, dt, parameters):
    v, s = states
    states[0] = v - s * dt
    states[1] = s + v * dt


def main():
    num_points = 5
    num_states = 2
    states = np.zeros((num_states, num_points))
    states[1, :] = np.linspace(0, 1, num_points)
    dt = 0.01
    t_bound = 20.0
    t0 = 0.0

    V_index = 0
    from time import perf_counter

    nT = int((t_bound - t0) / dt) - 1
    V = np.zeros((nT, num_points))
    t0 = perf_counter()
    beat.odesolver.solve(
        fun=simple_ode_forward_euler,
        t_bound=t_bound,
        states=states,
        V=V,
        V_index=V_index,
        dt=dt,
        parameters=None,
    )

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for i in range(num_points):
        ax.plot(V[:, i])
    fig.savefig("simple.png")


if __name__ == "__main__":
    main()
