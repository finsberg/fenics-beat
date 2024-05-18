# # Single-cell APD restitution
import matplotlib.pyplot as plt
import ap_features as apf
from pathlib import Path
import numpy as np
import gotranx
from beat.single_cell import get_steady_state, compute_hash

# import ap_features

here = Path(__file__).parent


model_path = Path("ORdmm_Land.py")
if not model_path.is_file():
    here = Path(__file__).parent
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

import ORdmm_Land

model = ORdmm_Land.__dict__


prebeats = 20
BCL = 1000
nbeats = 5
CI0 = 50
CI1 = BCL
CIinc = 25


outdir = here / "apd_restitution"
outdir.mkdir(exist_ok=True, parents=True)

track_indices = [model["state_index"]("v"), model["state_index"]("cai")]
save_every_ms = 2.0
dt = 0.01


# Pre-beats
y = get_steady_state(
    fun=model["forward_generalized_rush_larsen"],
    init_states=model["init_state_values"](),
    parameters=model["init_parameter_values"](celltype=0),
    outdir=outdir / "prebeats",
    track_indices=track_indices,
    nbeats=prebeats,
    save_every_ms=save_every_ms,
    BCL=BCL,
    dt=dt,
)

CIs = np.arange(CI1, CI0 - CIinc, -CIinc)
APDs = []
# APD restitution
for i, CI in enumerate(CIs):
    print("\nCI:", CI)
    out = outdir / f"CI_{CI}"
    y = get_steady_state(
        fun=model["forward_generalized_rush_larsen"],
        init_states=y,
        parameters=model["init_parameter_values"](celltype=0),
        outdir=out,
        track_indices=track_indices,
        nbeats=nbeats,
        save_every_ms=save_every_ms,
        BCL=CI,
        dt=dt,
    )
    traced_values_file = list(out.glob("tracked_values_*.npy"))[0]
    N = int(np.ceil(CI / save_every_ms))
    V = np.load(traced_values_file)[-N:, 0]
    t = np.linspace(0, CI, N)

    APD = apf.apd(factor=90, V=V, t=t)

    APDs.append(APD)
    print("APD:", APD)

fig, ax = plt.subplots()
ax.plot(CIs, APDs, "o-")
ax.set_xlabel("DI (ms)")
ax.set_ylabel("APD90 (ms)")
fig.savefig(outdir / "apd_restitution.png")
