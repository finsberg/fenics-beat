import logging
from typing import NamedTuple
import dolfin
import pint

try:
    import ufl_legacy as ufl
except ImportError:
    import ufl


from .units import ureg, to_quantity

logger = logging.getLogger(__name__)


def get_dimesion(u):
    # TODO : Check argument
    try:
        dim = ufl.domain.find_geometric_dimension(u)

    except ufl.UFLException as ex:
        try:
            dim = len(u)
        except Exception as ex2:
            logger.warning(ex)
            logger.warning(ex2)
            # Assume dimension is 3
            logger.warning("Assume dimension is 3")
            dim = 3

    return dim


def default_conductivities(name="Niederer") -> dict[str, pint.Quantity]:
    if name == "Niederer":
        return {
            "g_il": 0.17 * ureg("S/m"),
            "g_it": 0.019 * ureg("S/m"),
            "g_el": 0.62 * ureg("S/m"),
            "g_et": 0.24 * ureg("S/m"),
        }
    elif name == "Bishop":
        return {
            "g_il": 0.34 * ureg("S/m"),
            "g_it": 0.060 * ureg("S/m"),
            "g_el": 0.12 * ureg("S/m"),
            "g_et": 0.08 * ureg("S/m"),
        }
    elif name == "Potse":
        return {
            "g_il": 3.0 * ureg("S/m"),
            "g_it": 0.3 * ureg("S/m"),
            "g_el": 3.0 * ureg("S/m"),
            "g_et": 1.2 * ureg("S/m"),
        }
    else:
        raise ValueError(f"Unknown conductivity tensor {name}")


class Conductivities(NamedTuple):
    s_l: float
    s_t: float


def get_harmonic_mean_conductivity(
    chi: float,
    g_il: float = 0.17,
    g_it: float = 0.019,
    g_el: float = 0.62,
    g_et: float = 0.24,
) -> tuple[float, float]:
    # Conductivities as defined by page 4339 of Niederer benchmark
    sigma_il = to_quantity(g_il, "S/m")
    sigma_it = to_quantity(g_it, "S/m")
    sigma_el = to_quantity(g_el, "S/m")
    sigma_et = to_quantity(g_et, "S/m")

    # Compute monodomain approximation by taking harmonic mean in each
    # direction of intracellular and extracellular part
    def harmonic_mean(a, b):
        return a * b / (a + b)

    sigma_l = harmonic_mean(sigma_il, sigma_el)
    sigma_t = harmonic_mean(sigma_it, sigma_et)

    # Scale conducitivites by 1/(chi)
    s_l = (sigma_l / chi).to("uA/mV").magnitude
    s_t = (sigma_t / chi).to("uA/mV").magnitude
    return Conductivities(s_l, s_t)


def define_conductivity_tensor(
    chi: float,
    f0: dolfin.Constant | dolfin.Function,
    g_il: float = 0.17,
    g_it: float = 0.019,
    g_el: float = 0.62,
    g_et: float = 0.24,
):
    s_l, s_t = get_harmonic_mean_conductivity(chi, g_il, g_it, g_el, g_et)
    dim = get_dimesion(f0)

    # Define conductivity tensor
    return s_l * ufl.outer(f0, f0) + s_t * (ufl.Identity(dim) - ufl.outer(f0, f0))
