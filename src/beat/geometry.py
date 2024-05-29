from typing import NamedTuple
import numpy as np
import dolfin


class Geometry(NamedTuple):
    mesh: dolfin.Mesh
    ffun: dolfin.MeshFunction | None = None
    markers: dict[str, tuple[int, int]] | None = None
    f0: dolfin.Constant | dolfin.Function | None = None
    s0: dolfin.Constant | dolfin.Function | None = None
    n0: dolfin.Constant | dolfin.Function | None = None


def get_2D_slab_microstructure(
    transverse: bool = False,
) -> tuple[dolfin.Constant, dolfin.Constant]:
    if transverse:
        f0 = dolfin.Constant((0.0, 1.0))
        s0 = dolfin.Constant((1.0, 0.0))
    else:
        f0 = dolfin.Constant((1.0, 0.0))
        s0 = dolfin.Constant((0.0, 1.0))

    return f0, s0


def get_3D_slab_microstructure(
    transverse: bool = False,
) -> tuple[dolfin.Constant, dolfin.Constant, dolfin.Constant]:
    if transverse:
        f0 = dolfin.Constant((0.0, 0.0, 1.0))
        s0 = dolfin.Constant((1.0, 0.0, 0.0))
        n0 = dolfin.Constant((0.0, 1.0, 0.0))
    else:
        f0 = dolfin.Constant((1.0, 0.0, 0.0))
        s0 = dolfin.Constant((0.0, 1.0, 0.0))
        n0 = dolfin.Constant((0.0, 0.0, 1.0))

    return f0, s0, n0


def get_2D_slab_mesh(dx, Lx, Ly):
    return dolfin.RectangleMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0),
        dolfin.Point(Lx, Ly),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
    )


def get_3D_slab_mesh(dx, Lx, Ly, Lz=0.0, dim=2):
    return dolfin.BoxMesh(
        dolfin.MPI.comm_world,
        dolfin.Point(0.0, 0.0, 0.0),
        dolfin.Point(Lx, Ly, Lz),
        int(np.rint((Lx / dx))),
        int(np.rint((Ly / dx))),
        int(np.rint((Lz / dx))),
    )


def get_3D_slab_geometry(Lx, Ly, Lz, dx, transverse=False) -> Geometry:
    mesh = get_3D_slab_mesh(dx, Lx, Ly, Lz)
    f0, s0, n0 = get_3D_slab_microstructure(transverse)
    return Geometry(mesh=mesh, f0=f0, s0=s0, n0=n0)


def get_2D_slab_geometry(Lx, Ly, dx, transverse=False) -> Geometry:
    mesh = get_2D_slab_mesh(dx, Lx, Ly)
    f0, s0 = get_2D_slab_microstructure(transverse)
    return Geometry(mesh=mesh, f0=f0, s0=s0)
