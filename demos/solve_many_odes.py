# # Solve many ODEs
# In this example we solve how to use the ODEsolver to solve many ODEs at the same time.
# This might be useful if you need to run a grid search to find
# parameters that best fit your data.

import numpy as np
import beat
import matplotlib.pyplot as plt
from time import perf_counter


model = beat.cellmodels.beeler_reuter_1977
# model = beat.cellmodels.tentusscher_panfilov_2006.epi

init_states = model.init_state_values()
parameters = model.init_parameter_values()
parameters[model.parameter_indices("IstimAmplitude")] = 1.0
# parameters = model.init_parameters()
# parameters["IstimAmplitude"] = 1.0
num_points = 100
num_states = len(init_states)
states = np.zeros((num_states, num_points))
states.T[:] = init_states
dt = 0.1
t_bound = 1000.0
t0 = 0.0
# np.random.seed(1)

# index = model.parameter_indices("g_s")
# g_s = np.zeros((num_points))
# g_s[:] = parameters[index] * np.ones(num_points) + 0.0001 * np.random.random(
#     size=num_points
# )
# Maybe add support for this later
# parameters["g_s"] = gs
# model.parameters["a"] = np.zeros((num_points, 1))
# model.parameters["a"][:, 0] = 0.13 * np.ones(num_points) + 0.1 * np.random.random(
# size=num_points
# )

V_index = model.state_indices("V")

nT = int((t_bound - t0) / dt) - 1
V = np.zeros((nT, num_points))
t0 = perf_counter()
# values = np.zeros_like(ic)
beat.odesolver.solve(
    fun=model.forward_generalized_rush_larsen,
    t_bound=t_bound,
    states=states,
    V=V,
    V_index=V_index,
    dt=dt,
    parameters=parameters,
    # extra={"g_s": g_s},
)
el = perf_counter() - t0
print(f"Elapsed time: {el} seconds")

fig, ax = plt.subplots()
for i in range(10):
    ax.plot(V[:, i])
fig.savefig("cell.png")
