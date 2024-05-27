from typing import NamedTuple

import dolfin


class Geometry(NamedTuple):
    mesh: dolfin.Mesh
    ffun: dolfin.MeshFunction | None = None
    markers: dict[str, tuple[int, int]] | None = None
    f0: dolfin.Constant | dolfin.Function | None = None
    s0: dolfin.Constant | dolfin.Function | None = None
    n0: dolfin.Constant | dolfin.Function | None = None
