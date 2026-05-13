"""Jetsimpy interface module

Provides thin wrappers around the jetsimpy scientific core, mirroring the API of
``orphans.grb_interface`` but using the jetsimpy implementation.
"""

# Standard library imports
from copy import deepcopy

# Third‑party imports – may be unavailable in the CI environment.
import numpy as np
import jetsimpy
from astropy.cosmology import Planck18 as cosmo


# Local imports – these modules already exist in the package.
# pylint: disable=invalid-name,too-many-arguments,too-many-positional-arguments,unused-import
from .tools import flux_to_mag, get_wl_and_nu_band
# pylint: enable=invalid-name,too-many-arguments,too-many-positional-arguments,unused-import
from .jetsimpy_configs import GRB_BASE_PARAMS


def make_jet_light_curve(
    E0: float = 1.0e53,
    thetaObs: float = 0.05,
    thetaCore: float = 0.1,
    lf: float = 300,
    freq: float = 5.0e14,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute a GRB light curve using jetsimpy.

    The signature mirrors ``orphans.grb_interface.make_grb_light_curve`` so that
    callers can switch the emission core via a simple import change.

    Parameters
    ----------
    E0: isotropic‑equivalent energy (erg)
    thetaObs: observer angle (rad)
    thetaCore: jet half‑opening angle (rad)
    lf: Lorentz factor
    freq: observing frequency (Hz)

    Returns
    -------
    nu, t_sec, Fnu_Jy
        ``nu`` – frequency array (Hz)
        ``t_sec`` – time array (seconds)
        ``Fnu_Jy`` – flux density in Jansky (converted from jetsimpy's mJy)
    """
    # Build a parameter dictionary compatible with the jetsimpy API.
    P = deepcopy(GRB_BASE_PARAMS)
    P["Eiso"] = E0
    P["theta_v"] = thetaObs
    P["theta_c"] = thetaCore
    P["lf"] = lf

    # Time grid (seconds) – same geometric spacing as the afterglowpy version.
    t_sec = np.geomspace(1.0e3, 1.0e7, 100)

    # Frequency array – constant across the curve.
    nu = np.empty(t_sec.shape)
    nu[:] = freq

    # jetsimpy requires a Jet object; we follow the masson implementation.
    jet = jetsimpy.Jet(
        jetsimpy.PowerLaw(P["theta_c"], P["Eiso"], lf0=P["lf"]),
        P["A"],
        P["n0"],
        spread=False,
        grid=jetsimpy.ForwardJetRes(P["theta_c"], 129),
    )

    # Flux density returned by jetsimpy is in mJy.
    fnu = jet.FluxDensity(t_sec, nu, P)
    Fnu_Jy = fnu * 1.0e-3
    return nu, t_sec, Fnu_Jy


def make_jet_spectrum(
    E0: float = 1.0e53,
    z: float = 1.0,
    n0: float = 1.0,
    thetaObs: float = 0.05,
    thetaCore: float = 0.1,
    lf: float = 300,
    tday: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Compute a GRB spectral energy distribution using jetsimpy.

    Mirrors ``orphans.grb_interface.make_grb_spectrum`` but calls the jetsimpy
    backend.
    """
    P = deepcopy(GRB_BASE_PARAMS)
    P["Eiso"] = E0
    P["theta_v"] = thetaObs
    P["theta_c"] = thetaCore
    P["lf"] = lf
    P["n0"] = n0
    P["z"] = z
    # Distance in cm – convert from Mpc using astropy.
    P["d"] = cosmo.luminosity_distance(z).value * 3.086e24

    # Wavelength band.
    wl_full_band, freq_full_band = get_wl_and_nu_band()

    jet = jetsimpy.Jet(
        jetsimpy.PowerLaw(P["theta_c"], P["Eiso"], lf0=P["lf"]),
        P["A"],
        P["n0"],
        spread=False,
        grid=jetsimpy.ForwardJetRes(P["theta_c"], 129),
    )

    # Convert days to seconds for the jetsimpy call.
    t_sec = tday * 24 * 3600
    fnu = jet.FluxDensity(t_sec, freq_full_band, P)
    Fnu_Jy = fnu * 1.0e-3
    return wl_full_band, freq_full_band, tday, Fnu_Jy
