"""orphans is a library of functions to simulate and analyze orphan afterglow with Rubin LSST"""

# Export core symbols
from .grb_interface import make_grb_light_curve, make_grb_spectrum, dump_wl_Fnu_spectrum
from .tools import flux_to_mag, get_wl_and_nu_band
from .grb_configs import GRB_BASE_PARAMS
from .jetsimpy_interface import make_jet_light_curve, make_jet_spectrum

__all__ = [
    "make_grb_light_curve",
    "make_grb_spectrum",
    "dump_wl_Fnu_spectrum",
    "flux_to_mag",
    "get_wl_and_nu_band",
    "GRB_BASE_PARAMS",
    "make_jet_light_curve",
    "make_jet_spectrum",
]
