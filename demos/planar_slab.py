# # Planar slab
# This demo simulates a planar slab of cardiac tissue with a stimulus applied
# at the left end. The stimulus is applied for 2 ms and then turned off.
# ECGs are computed at 6 points on the surface of the slab.
#

import dolfin
import numpy as np
from collections import defaultdict
from pathlib import Path
import numpy.typing as npt

import matplotlib.pyplot as plt
from tqdm import tqdm

import beat

import beat.cellmodels.tentusscher_panfilov_2006 as model

# import beat.cellmodels.torord_dyn_chloride as model

model_name = model.__name__.split(".")[-1]


def setup_initial_conditions(func) -> npt.NDArray:
    ic = {
        "V": -85.23,  # mV
        "Xr1": 0.00621,
        "Xr2": 0.4712,
        "Xs": 0.0095,
        "m": 0.00172,
        "h": 0.7444,
        "j": 0.7045,
        "d": 3.373e-05,
        "f": 0.7888,
        "f2": 0.9755,
        "fCass": 0.9953,
        "s": 0.999998,
        "r": 2.42e-08,
        "Ca_i": 0.000126,  # millimolar
        "R_prime": 0.9073,
        "Ca_SR": 3.64,  # millimolar
        "Ca_ss": 0.00036,  # millimolar
        "Na_i": 8.604,  # millimolar
        "K_i": 136.89,  # millimolar
    }
    return func(**ic)


def setup_geometry(dx):
    Lx = 20.0  # mm
    Ly = 7.0  # mm
    Lz = 3.0  # mm

    mesh = dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
        int(np.rint((Lz / dx))),
    )

    return mesh


def define_stimulus(mesh, chi, C_m, time):
    S1_marker = 1
    L = 1.5
    S1_subdomain = dolfin.CompiledSubDomain(
        "x[0] <= L + DOLFIN_EPS",
        L=L,
    )
    S1_markers = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
    S1_subdomain.mark(S1_markers, S1_marker)

    A = 50000.0  # mu A/cm^3
    cm2mm = 10.0
    factor = 1.0 / (chi * C_m)  # NB: cbcbeat convention
    amplitude = factor * A * (1.0 / cm2mm) ** 3  # mV/ms

    duration = 2
    I_s = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=0.0,
        duration=duration,
        amplitude=amplitude,
        degree=0,
    )
    dx = dolfin.Measure("dx", domain=mesh, subdomain_data=S1_markers)(S1_marker)
    return beat.base_model.Stimulus(dz=dx, expr=I_s)


def define_conductivity_tensor(chi, C_m):
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = 0.17  # mS / mm
    sigma_it = 0.019  # mS / mm
    sigma_el = 0.62  # mS / mm
    sigma_et = 0.24  # mS / mm

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(C_m * chi)
    s_l = sigma_l / (C_m * chi)  # mm^2 / ms
    s_t = sigma_t / (C_m * chi)  # mm^2 / ms

    # Define conductivity tensor
    M = dolfin.as_tensor(((s_l, 0, 0), (0, s_t, 0), (0, 0, s_t)))

    return M


def load_timesteps_from_xdmf(xdmffile):
    import xml.etree.ElementTree as ET

    times = {}
    i = 0
    tree = ET.parse(xdmffile)
    for elem in tree.iter():
        if elem.tag == "Time":
            times[i] = float(elem.get("Value"))
            i += 1

    return times


def load_from_file(heart_mesh, xdmffile, key="v", stop_index=None):
    V = dolfin.FunctionSpace(heart_mesh, "Lagrange", 1)
    v = dolfin.Function(V)

    timesteps = load_timesteps_from_xdmf(xdmffile)
    with dolfin.XDMFFile(Path(xdmffile).as_posix()) as f:
        for i, t in tqdm(timesteps.items()):
            f.read_checkpoint(v, key, i)
            yield v.copy(deepcopy=True), t


def compute_ecg_recovery():
    datadir = Path("planar_slab")
    xdmffile = datadir / f"state_{model_name}.xdmf"
    mesh = dolfin.Mesh()

    with dolfin.XDMFFile((datadir / "mesh.xdmf").as_posix()) as xdmf:
        xdmf.read(mesh)

    vs = load_from_file(mesh, xdmffile, key="V")

    leads = dict(
        p1=(22.0, 8.0, 4.0),
        p2=(22.0, 0.0, 4.0),
        p3=(22.0, 8.0, 0.0),
        p4=(22.0, -1.0, -1.0),
        p5=(22.0, 3.5, 4.0),
        p6=(22.0, 3.5, -1.0),
        ground=(-10.0, -10.0, -10.0),
    )

    fname = datadir / f"extracellular_potential_{model_name}.npy"
    if not fname.is_file():
        phie = defaultdict(list)
        time = []
        for v, t in vs:
            time.append(t)
            for name, point in leads.items():
                phie[name].append(
                    beat.ecg.ecg_recovery(v=v, mesh=mesh, sigma_b=1.0, point=point)
                )
        np.save(fname, {"phie": phie, "time": time})

    phie_time = np.load(fname, allow_pickle=True).item()
    phie = phie_time["phie"]
    time = phie_time["time"]

    fig, ax = plt.subplots(2, 4, sharex=True, figsize=(12, 8))
    for i, (name, values) in enumerate(phie.items()):
        axi = ax.flatten()[i]
        axi.plot(time, values)
        axi.set_title(name)
    fig.savefig(datadir / f"extracellular_potential_{model_name}.png")

    # ecg = beat.ecg.Leads12(**{k: np.array(v) for k, v in phie.items()})
    fig, ax = plt.subplots()
    for i, name in enumerate(["p1", "p2", "p3", "p4", "p5", "p6"]):
        y = np.subtract(phie[name], phie["ground"])
        ax.plot(time, y, label=name)

    fig.savefig(datadir / f"ecg_{model_name}.png")


def main(datadir=Path("planar_slab")):
    Lx = 20.0  # mm
    Ly = 7.0  # mm
    Lz = 3.0  # mm
    dx = 0.5

    mesh = dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
        int(np.rint((Lz / dx))),
    )
    with dolfin.XDMFFile((datadir / "mesh.xdmf").as_posix()) as xdmf:
        xdmf.write(mesh)

    mfun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    mfun.set_all(0)
    dolfin.CompiledSubDomain("near(x[0], 0)").mark(mfun, 1)
    dolfin.CompiledSubDomain(f"near(x[0], {Lx})").mark(mfun, 2)

    V = dolfin.FunctionSpace(mesh, "Lagrange", 1)
    markers = dolfin.Function(V)
    arr = markers.vector().get_local().copy()
    v2d = dolfin.vertex_to_dof_map(V)

    for marker in [0, 1, 2]:
        for facet in mfun.where_equal(marker):
            f = dolfin.Facet(mesh, facet)
            arr[v2d[f.entities(0)]] = marker
    markers.vector().set_local(arr)

    with dolfin.XDMFFile((datadir / "markers.xdmf").as_posix()) as xdmf:
        xdmf.write(markers)

    # Use some custom ic for tentusscher
    setup_ic = (
        setup_initial_conditions
        if model_name == "tentusscher_panfilov_2006"
        else lambda f: f()
    )

    init_states = {
        0: setup_ic(model.mid.init_state_values),
        1: setup_ic(model.endo.init_state_values),
        2: setup_ic(model.epi.init_state_values),
    }
    parameters = {
        0: model.mid.init_parameter_values(),
        1: model.endo.init_parameter_values(),
        2: model.epi.init_parameter_values(),
    }
    fun = {
        0: model.mid.forward_generalized_rush_larsen,
        1: model.endo.forward_generalized_rush_larsen,
        2: model.epi.forward_generalized_rush_larsen,
    }
    if model_name == "tentusscher_panfilov_2006":
        v_index = {
            0: model.mid.state_indices("V"),
            1: model.endo.state_indices("V"),
            2: model.epi.state_indices("V"),
        }
        parameters[0][model.mid.parameter_indices("stim_amplitude")] = 0.0
        parameters[1][model.endo.parameter_indices("stim_amplitude")] = 0.0
        parameters[2][model.epi.parameter_indices("stim_amplitude")] = 0.0

    elif model_name == "torord_dyn_chloride":
        v_index = {
            0: model.mid.state_indices("v"),
            1: model.endo.state_indices("v"),
            2: model.epi.state_indices("v"),
        }
        parameters[0][model.mid.parameter_indices("i_Stim_Amplitude")] = 0.0
        parameters[1][model.endo.parameter_indices("i_Stim_Amplitude")] = 0.0
        parameters[2][model.epi.parameter_indices("i_Stim_Amplitude")] = 0.0
    # Surface to volume ratio
    chi = 140.0  # mm^{-1}
    # Membrane capacitance
    C_m = 0.01  # mu F / mm^2

    time = dolfin.Constant(0.0)
    I_s = define_stimulus(
        mesh=mesh,
        chi=chi,
        C_m=C_m,
        time=time,
    )

    M = define_conductivity_tensor(chi, C_m)

    params = {"preconditioner": "sor", "use_custom_preconditioner": False}
    pde = beat.MonodomainModel(time=time, mesh=mesh, M=M, I_s=I_s, params=params)

    ode = beat.odesolver.DolfinMultiODESolver(
        pde.state,
        markers=markers,
        num_states={i: len(s) for i, s in init_states.items()},
        fun=fun,
        init_states=init_states,
        parameters=parameters,
        v_index=v_index,
    )

    T = 1
    # Change to 500 to simulate the full cardiac cycle
    # T = 500
    t = 0.0
    dt = 0.05
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

    fname = (datadir / f"state_{model_name}.xdmf").as_posix()
    i = 0
    while t < T + 1e-12:
        if i % 20 == 0:
            v = solver.pde.state.vector().get_local()
            print(f"Solve for {t=:.2f}, {v.max() =}, {v.min() = }")
            with dolfin.XDMFFile(dolfin.MPI.comm_world, fname) as xdmf:
                xdmf.write_checkpoint(
                    solver.pde.state,
                    "V",
                    float(t),
                    dolfin.XDMFFile.Encoding.HDF5,
                    True,
                )
        solver.step((t, t + dt))
        i += 1
        t += dt


if __name__ == "__main__":
    # main()
    compute_ecg_recovery()
