# # Conduction velocity and ECG for slabs
# In this demo we will show how to compute conduction velocity and ECG for a Slab geometry.
#

import argparse
from typing import Sequence, NamedTuple
from pathlib import Path
import beat.single_cell
import dolfin
import matplotlib.pyplot as plt
import numpy as np
import pint
import gotranx
import beat


try:
    import ufl_legacy as ufl
except ImportError:
    import ufl

ureg = pint.UnitRegistry()


class Geometry(NamedTuple):
    mesh: dolfin.Mesh
    ffun: dolfin.MeshFunction
    markers: dict[str, tuple[int, int]]
    f0: dolfin.Constant
    s0: dolfin.Constant
    n0: dolfin.Constant


def define_stimulus(
    mesh, chi, time, ffun, markers, mesh_unit, duration=2.0, stim_amp=500.0
):
    duration = 2.0  # ms
    A = stim_amp * ureg("uA/cm**2")

    amplitude = (A / chi).to(f"uA/{mesh_unit}").magnitude

    I_s = dolfin.Expression(
        "time >= start ? (time <= (duration + start) ? amplitude : 0.0) : 0.0",
        time=time,
        start=0.0,
        duration=duration,
        amplitude=amplitude,
        degree=0,
    )

    subdomain_data = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
    subdomain_data.set_all(0)
    marker = 1
    subdomain_data.array()[ffun.array() == markers["ENDO"][0]] = 1

    ds = dolfin.Measure("ds", domain=mesh, subdomain_data=subdomain_data)(marker)
    return beat.base_model.Stimulus(dz=ds, expr=I_s)


def default_conductivities(name="Niederer") -> dict[str, float]:
    if name == "Niederer":
        return {
            "g_il": 0.17,
            "g_it": 0.019,
            "g_el": 0.62,
            "g_et": 0.24,
        }
    elif name == "Bishop":
        return {
            "g_il": 0.34,
            "g_it": 0.060,
            "g_el": 0.12,
            "g_et": 0.08,
        }
    elif name == "Potse":
        return {
            "g_il": 3.0,
            "g_it": 0.3,
            "g_el": 3.0,
            "g_et": 1.2,
        }
    else:
        raise ValueError(f"Unknown conductivity tensor {name}")


def define_conductivity_tensor(
    chi,
    f0,
    s0,
    n0,
    g_il=0.17,
    g_it=0.019,
    g_el=0.62,
    g_et=0.24,
    dim: int = 2,
):
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = g_il * ureg("S/m")
    sigma_it = g_it * ureg("S/m")
    sigma_el = g_el * ureg("S/m")
    sigma_et = g_et * ureg("S/m")

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(chi)
    s_l = (sigma_l / chi).to("uA/mV").magnitude
    s_t = (sigma_t / chi).to("uA/mV").magnitude

    # Define conductivity tensor
    if dim == 2:
        A = dolfin.as_matrix(
            [
                [f0[0], s0[0]],
                [f0[1], s0[1]],
            ],
        )

        M_star = ufl.diag(dolfin.as_vector([s_l, s_t]))
    else:
        A = dolfin.as_matrix(
            [
                [f0[0], s0[0], n0[0]],
                [f0[1], s0[1], n0[1]],
                [f0[2], s0[2], n0[2]],
            ],
        )

        M_star = ufl.diag(dolfin.as_vector([s_l, s_t, s_t]))
    M = A * M_star * A.T

    return M


def setup_geometry(dx, Lx, Ly, Lz=0.0, dim=2):
    if dim == 2:
        mesh = dolfin.RectangleMesh(
            dolfin.MPI.comm_world,
            dolfin.Point(0.0, 0.0),
            dolfin.Point(Lx, Ly),
            int(np.rint((Lx / dx))),
            int(np.rint((Ly / dx))),
        )

    else:
        mesh = dolfin.BoxMesh(
            dolfin.MPI.comm_world,
            dolfin.Point(0.0, 0.0, 0.0),
            dolfin.Point(Lx, Ly, Lz),
            int(np.rint((Lx / dx))),
            int(np.rint((Ly / dx))),
            int(np.rint((Lz / dx))),
        )
    return mesh


def compute_conduction_velocity(
    mesh, results_folder, L, dx, mesh_unit: str, threshold=0.0
):
    x0 = L * 0.25
    x1 = L * 0.75
    if mesh.geometry().dim() == 2:
        p1 = (x0, dx * 0.5)
        p2 = (x1, dx * 0.5)
    else:
        p1 = (x0, dx * 0.5, dx * 0.5)  # type: ignore
        p2 = (x1, dx * 0.5, dx * 0.5)  # type: ignore

    V = dolfin.FunctionSpace(mesh, "CG", 1)
    v = dolfin.Function(V)

    fname = (results_folder / "V.xdmf").as_posix()
    times = load_timesteps_from_xdmf(fname)

    t1 = np.inf
    t2 = np.inf

    for i, t in times.items():
        with dolfin.XDMFFile(mesh.mpi_comm(), fname) as xdmf:
            xdmf.read_checkpoint(v, "V", i)
            print(f"Read {t=:.2f}, {v(p1) =}, {v(p2) = }")

        if v(p1) > threshold:
            t1 = min(t, t1)
        if v(p2) > threshold:
            t2 = min(t, t2)
            break

    cv = (x1 - x0) / (t2 - t1) * ureg(f"{mesh_unit}/ms")
    msg = (
        f"Conduction velocity = {cv.magnitude:.3f} mm/ms or "  #
        f" {cv.to('m/s').magnitude:.3f} m/s or "  #
        f" {cv.to('cm/s').magnitude:.3f} cm/s"  #
    )
    print(msg)


def compute_ecg_recovery(mesh, results_folder, L, dx):
    x0 = L * 2.0
    if mesh.geometry().dim() == 2:
        p = (x0, dx * 0.5)
    else:
        p = (x0, dx * 0.5, dx * 0.5)

    V = dolfin.FunctionSpace(mesh, "CG", 1)
    v = dolfin.Function(V)

    fname = (results_folder / "V.xdmf").as_posix()
    times = load_timesteps_from_xdmf(fname)

    phie = []
    for i, t in times.items():
        with dolfin.XDMFFile(mesh.mpi_comm(), fname) as xdmf:
            xdmf.read_checkpoint(v, "V", i)
            print(f"Read {t=:.2f}")
            phie.append(beat.ecg.ecg_recovery(v=v, mesh=mesh, sigma_b=1.0, point=p))

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.plot(times.values(), phie)
    ax.set_title("ECG recovery")
    fig.savefig(results_folder / "ecg_recovery.png")


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


def delete_old_file(file: str):
    Path(file).unlink(missing_ok=True)
    Path(file).with_suffix(".h5").unlink(missing_ok=True)


def run_model(
    data: Geometry,
    markers: dolfin.Function,
    resultsdir: Path,
    end_time: float = 15.0,
    dt: float = 0.05,
    mesh_unit="cm",
    save_freq: int = 20,
    stim_amp: float = 5000.0,
    g_il=0.16069,
    g_el=0.625,
    g_it=0.04258,
    g_et=0.236,
    **kwargs,
):
    print("Running model")
    # Load the model
    model_path = Path("ORdmm_Land.py")
    if not model_path.is_file():
        print("Generate code for cell model")
        here = Path.cwd()
        ode = gotranx.load_ode(here / "ORdmm_Land.ode")
        code = gotranx.cli.gotran2py.get_code(
            ode, scheme=[gotranx.schemes.Scheme.forward_generalized_rush_larsen]
        )
        model_path.write_text(code)

    import ORdmm_Land

    model = ORdmm_Land.__dict__

    # Surface to volume ratio
    chi = 1400.0 * ureg("cm**-1")

    # Membrane capacitance
    C_m = 1.0 * ureg("uF/cm**2")

    with dolfin.XDMFFile((resultsdir / "markers.xdmf").as_posix()) as xdmf:
        xdmf.write(markers)

    print("Get steady states")
    nbeats = 2  # Should be set to at least 200
    init_states = {
        0: beat.single_cell.get_steady_state(
            fun=model["forward_generalized_rush_larsen"],
            init_states=model["init_state_values"](),
            parameters=model["init_parameter_values"](celltype=2),
            outdir=resultsdir / "mid",
            BCL=1000,
            nbeats=nbeats,
            track_indices=[model["state_index"]("v"), model["state_index"]("cai")],
            dt=0.05,
        ),
        1: beat.single_cell.get_steady_state(
            fun=model["forward_generalized_rush_larsen"],
            init_states=model["init_state_values"](),
            parameters=model["init_parameter_values"](celltype=0),
            outdir=resultsdir / "endo",
            BCL=1000,
            nbeats=nbeats,
            track_indices=[
                model["state_index"]("v"),
                model["state_index"]("cai"),
                model["state_index"]("nai"),
            ],
            dt=0.05,
        ),
        2: beat.single_cell.get_steady_state(
            fun=model["forward_generalized_rush_larsen"],
            init_states=model["init_state_values"](),
            parameters=model["init_parameter_values"](celltype=1),
            outdir=resultsdir / "epi",
            BCL=1000,
            nbeats=nbeats,
            track_indices=[model["state_index"]("v"), model["state_index"]("cai")],
            dt=0.05,
        ),
    }
    # endo = 0, epi = 1, M = 2
    parameters = {
        0: model["init_parameter_values"](amp=0.0, celltype=2),
        1: model["init_parameter_values"](amp=0.0, celltype=0),
        2: model["init_parameter_values"](amp=0.0, celltype=1),
    }
    fun = {
        0: model["forward_generalized_rush_larsen"],
        1: model["forward_generalized_rush_larsen"],
        2: model["forward_generalized_rush_larsen"],
    }
    v_index = {
        0: model["state_index"]("v"),
        1: model["state_index"]("v"),
        2: model["state_index"]("v"),
    }

    time = dolfin.Constant(0.0)
    I_s = define_stimulus(
        mesh=data.mesh,
        chi=chi,
        time=time,
        ffun=data.ffun,
        markers=data.markers,
        stim_amp=stim_amp,
        mesh_unit=mesh_unit,
    )

    M = define_conductivity_tensor(
        chi,
        f0=data.f0,
        s0=data.s0,
        n0=data.n0,
        g_il=g_il,
        g_it=g_it,
        g_el=g_el,
        g_et=g_et,
        dim=data.mesh.geometry().dim(),
    )

    params = {"preconditioner": "sor", "use_custom_preconditioner": False}
    pde = beat.MonodomainModel(
        time=time,
        mesh=data.mesh,
        M=M,
        I_s=I_s,
        params=params,
        C_m=C_m.to(f"uF/{mesh_unit}**2").magnitude,
    )

    V_ode = dolfin.FunctionSpace(data.mesh, "Lagrange", 1)
    ode = beat.odesolver.DolfinMultiODESolver(
        v_ode=dolfin.Function(V_ode),
        v_pde=pde.state,
        markers=markers,
        num_states={i: len(s) for i, s in init_states.items()},
        fun=fun,
        init_states=init_states,
        parameters=parameters,
        v_index=v_index,
    )
    t = 0.0
    solver = beat.MonodomainSplittingSolver(pde=pde, ode=ode)

    fname = (resultsdir / "V.xdmf").as_posix()
    delete_old_file(fname)

    def save(t):
        v = solver.pde.state.vector().get_local()
        print(f"Solve for {t=:.2f}, {v.max() =}, {v.min() = }")
        with dolfin.XDMFFile(data.mesh.mpi_comm(), fname) as xdmf:
            xdmf.parameters["functions_share_mesh"] = True
            xdmf.parameters["rewrite_function_mesh"] = False
            xdmf.write_checkpoint(
                solver.pde.state,
                "V",
                float(t),
                dolfin.XDMFFile.Encoding.HDF5,
                True,
            )

    i = 0
    while t < end_time + 1e-12:
        # Make sure to save at the same time steps that is used by Ambit

        if i % save_freq == 0:
            save(t)

        solver.step((t, t + dt))
        i += 1
        t += dt


def get_microstructure(
    dim: int, transverse: bool = False
) -> tuple[dolfin.Constant, dolfin.Constant, dolfin.Constant]:
    if dim == 2:
        if transverse:
            f0 = dolfin.Constant((0.0, 1.0))
            s0 = dolfin.Constant((1.0, 0.0))
        else:
            f0 = dolfin.Constant((1.0, 0.0))
            s0 = dolfin.Constant((0.0, 1.0))

        n0 = dolfin.Constant((0.0, 0.0))

    else:
        if transverse:
            f0 = dolfin.Constant((0.0, 0.0, 1.0))
            s0 = dolfin.Constant((1.0, 0.0, 0.0))
            n0 = dolfin.Constant((0.0, 1.0, 0.0))
        else:
            f0 = dolfin.Constant((1.0, 0.0, 0.0))
            s0 = dolfin.Constant((0.0, 1.0, 0.0))
            n0 = dolfin.Constant((0.0, 0.0, 1.0))

    return f0, s0, n0


def main() -> int:
    results_folder = Path("results-slab")
    save_every_ms = 1.0
    dimension = 2
    transverse = False
    end_time = 20.0
    dt = 0.05
    overwrite = False
    stim_amp = 5000.0
    mesh_unit = "cm"

    # Load mesh
    mesh_unit = mesh_unit

    dx = 0.05 * ureg("cm").to(mesh_unit).magnitude
    L = 1.0 * ureg("cm").to(mesh_unit).magnitude
    mesh = setup_geometry(Lx=L, Ly=dx, Lz=dx, dx=dx / 5, dim=dimension)
    if not results_folder.is_dir() or overwrite:
        ffun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim() - 1)
        ffun.set_all(0)
        stim_domain = dolfin.CompiledSubDomain("x[0] <= DOLFIN_EPS")
        marker = 1
        stim_domain.mark(ffun, marker)

        cfun = dolfin.MeshFunction("size_t", mesh, mesh.topology().dim())
        cfun.set_all(0)
        endo = dolfin.CompiledSubDomain("x[0] < L / 3", L=L)
        epi = dolfin.CompiledSubDomain("x[0] > 2 * L / 3", L=L)
        endo.mark(cfun, 1)
        epi.mark(cfun, 2)
        V = dolfin.FunctionSpace(mesh, "CG", 1)
        cfun_DG = dolfin.Function(dolfin.FunctionSpace(mesh, "DG", 0))
        cfun_DG.vector()[:] = cfun.array()
        cfun_func = dolfin.Function(V)
        cfun_func.interpolate(cfun_DG)

        g_il = 0.16069
        g_el = 0.625
        g_it = 0.04258
        g_et = 0.236

        f0, s0, n0 = get_microstructure(dimension, transverse)

        markers = {"ENDO": (marker, 2)}

        data = Geometry(
            mesh=mesh,
            ffun=ffun,
            markers=markers,
            f0=f0,
            s0=s0,
            n0=n0,
        )
        save_freq = round(save_every_ms / dt)
        run_model(
            data=data,
            markers=cfun_func,
            g_il=g_il,
            g_it=g_it,
            g_el=g_el,
            g_et=g_et,
            save_freq=save_freq,
            resultsdir=results_folder,
            end_time=end_time,
            dt=dt,
            stim_amp=stim_amp,
            mesh_unit=mesh_unit,
        )
    else:
        print(f"Results folder {results_folder} exists")

    compute_conduction_velocity(mesh, results_folder, L, dx, mesh_unit)
    compute_ecg_recovery(mesh, results_folder, L, dx)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
