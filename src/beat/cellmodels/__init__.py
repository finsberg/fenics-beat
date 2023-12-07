from . import beeler_reuter_1977
from . import fitzhughnagumo
from . import tentusscher_panfilov_2006
from . import torord_dyn_chloride


def all_cellmodels():
    return [
        fitzhughnagumo,
        beeler_reuter_1977,
        tentusscher_panfilov_2006.endo,
        tentusscher_panfilov_2006.epi,
        tentusscher_panfilov_2006.mid,
        torord_dyn_chloride.endo,
        torord_dyn_chloride.epi,
        torord_dyn_chloride.mid,
    ]


__all__ = [
    "fitzhughnagumo",
    "beeler_reuter_1977",
    "tentusscher_panfilov_2006",
    "torord_dyn_chloride",
]
