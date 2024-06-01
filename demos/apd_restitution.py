# # Single-cell APD restitution
# TBW

# +
import matplotlib.pyplot as plt
import ap_features as apf
from pathlib import Path
import numpy as np
import gotranx
from beat.single_cell import get_steady_state, compute_hash


here = Path.cwd()


# +
model_path = Path("tentusscher_panfilov_2006_epi_cell.py")
if not model_path.is_file():
    ode = gotranx.load_ode(
        here
        / ".."
        / "odes"
        / "tentusscher_panfilov_2006"
        / "tentusscher_panfilov_2006_epi_cell.ode"
    )
    code = gotranx.cli.gotran2py.get_code(
        ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
    )
    model_path.write_text(code)

import tentusscher_panfilov_2006_epi_cell

model = tentusscher_panfilov_2006_epi_cell.__dict__


# -


def run(outdir, parameters):

    # Pre-beats
    y = get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=model["init_state_values"](),
        parameters=parameters,
        outdir=outdir / "prebeats",
        track_indices=track_indices,
        nbeats=prebeats,
        save_every_ms=save_every_ms,
        BCL=BCL,
        dt=dt,
    )

    CIs = np.arange(CI1, CI0 - CIinc, -CIinc)
    APDs = []
    Vs = []
    stim_period_index = model["parameter_index"]("stim_period")
    # APD restitution
    for CI in CIs:
        parameters[stim_period_index] = CI
        print("\nCI:", CI)
        out = outdir / f"CI_{CI}"
        hash_input = compute_hash(
            fun=model["forward_generalized_rush_larsen"],
            init_states=y,
            parameters=parameters,
            nbeats=nbeats,
            BCL=CI,
            dt=dt,
        )
        y = get_steady_state(
            fun=model["forward_generalized_rush_larsen"],
            init_states=y,
            parameters=parameters,
            outdir=out,
            track_indices=track_indices,
            nbeats=nbeats,
            save_every_ms=save_every_ms,
            BCL=CI,
            dt=dt,
        )

        traced_values_file = out / f"tracked_values_{hash_input}.npy"
        N = int(np.ceil(CI / save_every_ms))

        V = np.load(traced_values_file)[-N:, 0]
        t = np.linspace(0, CI, N)
        Vs.append(V)

        APD = apf.apd(factor=90, V=V[-N:], t=t)

        APDs.append(APD)
        print("APD:", APD)

    return APDs, CIs, Vs


prebeats = 1  # Use 50 or more
BCL = 1000
nbeats = 1  # Use maybe 7 - 8 beats here
CI0 = 200  # Final interval - typically 200 ms
# How much to decrease the interval by each step
# Typical value here is 25
CIinc = 100
CI1 = BCL


outdir = here / "apd_restitution"
outdir.mkdir(exist_ok=True, parents=True)

track_indices = [model["state_index"]("V"), model["state_index"]("Ca_i")]
save_every_ms = 2.0
dt = 0.05


APDs_normal, CIs, Vs = run(
    outdir=outdir / "normal",
    parameters=model["init_parameter_values"](),
)
fig, ax = plt.subplots()
ax.plot(CIs, APDs_normal, "o-")
ax.set_xlabel("DI (ms)")
ax.set_ylabel("APD90 (ms)")
fig.savefig(outdir / "apd_restitution_normal.png")

step = 3
fig_v, ax_v = plt.subplots(1, 3, sharex=True, sharey=True)

lines = []
labels = []
for i, V in enumerate(Vs[::step]):
    j = i * step
    (l,) = ax_v[0].plot(np.linspace(0, CIs[j], len(V)), V)
    lines.append(l)
    labels.append(f"{CIs[j]} ms")
ax_v[0].set_title("Normal")

APDs_step1, CIs, Vs = run(
    outdir=outdir / "step1",
    parameters=model["init_parameter_values"](g_Kr=0.172, g_Ks=0.441),
)
fig, ax = plt.subplots()
ax.plot(CIs, APDs_step1, "o-")
ax.set_xlabel("DI (ms)")
ax.set_ylabel("APD90 (ms)")
fig.savefig(outdir / "apd_restitution_step1.png")

for i, V in enumerate(Vs[::step]):
    j = i * step
    ax_v[1].plot(np.linspace(0, CIs[j], len(V)), V)
ax_v[1].set_title("Step 1")

APDs_step2, CIs, Vs = run(
    outdir=outdir / "step2",
    parameters=model["init_parameter_values"](g_Kr=0.134, g_Ks=0.270),
)
fig, ax = plt.subplots()
ax.plot(CIs, APDs_step2, "o-")
ax.set_xlabel("DI (ms)")
ax.set_ylabel("APD90 (ms)")
fig.savefig(outdir / "apd_restitution_step2.png")

for i, V in enumerate(Vs[::step]):
    j = i * step
    ax_v[2].plot(np.linspace(0, CIs[j], len(V)), V)
ax_v[2].set_title("Step 2")

fig, ax = plt.subplots()
ax.plot(CIs, APDs_normal, "o-", label="normal")
ax.plot(CIs, APDs_step1, "o-", label="step1")
ax.plot(CIs, APDs_step2, "o-", label="step2")
ax.set_xlabel("DI (ms)")
ax.set_ylabel("APD90 (ms)")
ax.legend()
fig.savefig(outdir / "apd_restitution.png")


for axi in ax_v:
    axi.set_xlabel("Time (ms)")
    axi.grid()

ax_v[0].set_ylabel("V (mV)")
lgd = fig_v.legend(lines, labels, loc="center left", bbox_to_anchor=(1, 0.5))
fig_v.tight_layout()
fig_v.savefig(outdir / "V.png", bbox_extra_artists=(lgd,), bbox_inches="tight")
