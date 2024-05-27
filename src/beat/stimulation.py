from typing import NamedTuple
import dolfin


from .units import ureg


class Stimulus(NamedTuple):
    dz: dolfin.Measure
    expr: dolfin.Expression


def define_stimulus(
    mesh: dolfin.Mesh,
    chi: float,
    time: dolfin.Constant,
    subdomain_data: dolfin.MeshFunction,
    marker: int,
    mesh_unit: str = "cm",
    duration: float = 2.0,
    amplitude: float = 500.0,
    start: float = 0.0,
    PCL: float | dolfin.Constant = 1000.0,
):
    dim = subdomain_data.dim()

    if isinstance(amplitude, ureg.Quantity):
        A = amplitude
    else:
        if dim == 1:
            A = amplitude * ureg("uA / cm")
        elif dim == 2:
            A = amplitude * ureg("uA / cm**2")
        elif dim == 3:
            A = amplitude * ureg("uA / cm**3")

    amp = (A / chi).to(f"uA/{mesh_unit}**{dim - 1}").magnitude

    I_s = dolfin.Expression(
        "std::fmod(time,PCL) >= start "
        "? (std::fmod(time,PCL) <= (duration + start) ? amplitude : 0.0)"
        " : 0.0",
        time=time,
        start=start,
        duration=duration,
        amplitude=amp,
        degree=0,
        PCL=PCL,
    )

    if dim == 1:
        dz = dolfin.Measure("dP", domain=mesh, subdomain_data=subdomain_data)(marker)
    elif dim == 2:
        dz = dolfin.Measure("ds", domain=mesh, subdomain_data=subdomain_data)(marker)
    elif dim == 3:
        dz = dolfin.Measure("dx", domain=mesh, subdomain_data=subdomain_data)(marker)
    return Stimulus(dz=dz, expr=I_s)
