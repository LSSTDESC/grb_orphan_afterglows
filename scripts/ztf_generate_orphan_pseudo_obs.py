import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import afterglowpy as grb
import random
import pickle

from astropy.cosmology import Planck18 as cosmo
from tqdm import tqdm

import skysurvey
import sncosmo

from skysurvey.target.core import Transient
from skysurvey.tools import random_radec

import warnings
warnings.filterwarnings("ignore")


#=======================================================================================================
def get_grb_params():
    """ Get an orphan parameter configuration among orphan that were observable by Rubin LSST

    Returns
    -------
    grb_params : dictionary
        Model parameters used to compute GRB afterglow light curves
    """

    return df_orphans_ztf.iloc[i]['config']


def get_flux(grb_params):
    """ Compute the spectral flux density for ZTF frequency coverage.
    The flux is computed by `afterglowpy` (Ryan et al. 2020).

    Parameters
    ----------
    grb_params : dictionary
        Model parameters used to compute GRB afterglow light curves

    Returns
    -------
    phase : numpy.ndarray
        Phases in days.
    wave : numpy.ndarray
        Wavelengths in Angstroms.
    flux : numpy.ndarray
        Model spectral flux density in arbitrary units.
    """

    phase = np.linspace(0., 100., 100)  # time in days
    wave = np.linspace(3000, 9000, 100)  # wavelength in Angstrom
    nu = 299792458 / (wave * 1e-10)

    flux = []
    for p in phase:
        if p == 0.:
            flux.append([10 ** -28] * len(phase))
        else:
            t = p * np.ones(nu.shape) * grb.day2sec
            mJys = grb.fluxDensity(t, nu, **grb_params)
            Jys = 1e-3 * mJys

            # convert to erg/s/cm^2/A
            flux.append(Jys * 1e-23 * 2.99792458e18 / (wave ** 2))  # 2.99792458e18 is the light speed in A/s

    return phase, wave, flux


def get_orphan_model():
    """ """

    grb_params = get_grb_params()
    phase, wave, flux = get_flux(grb_params)

    orphan_source = sncosmo.TimeSeriesSource(phase, wave, np.array(flux), name='orphan')
    model = sncosmo.Model(orphan_source)

    return model


def save_ztf_light_curve(grb_params, transient_data, dset_data, index_max):
    # keep only the times, flux and zero-points 5 days before and 50 days after the beginning of the light curve
    tmin = transient_data['t0'] - 5.
    tmax = transient_data['t0'] + 50.

    time_obs = np.array(dset_data.loc[(dset_data['time'] >= tmin) & (dset_data['time'] <= tmax), 'time'])
    flux_obs = np.array(dset_data.loc[(dset_data['time'] >= tmin) & (dset_data['time'] <= tmax), 'flux'])
    flux_err_obs = np.array(dset_data.loc[(dset_data['time'] >= tmin) & (dset_data['time'] <= tmax), 'fluxerr'])
    zp_obs = np.array(dset_data.loc[(dset_data['time'] >= tmin) & (dset_data['time'] <= tmax), 'zp'])
    filts = np.array(dset_data.loc[(dset_data['time'] >= tmin) & (dset_data['time'] <= tmax), 'band'])

    # convert from observed flux to magnitude using the zero-point value of the observation
    mag_obs = -2.5 * np.log10(flux_obs) + zp_obs

    infos = {'config': grb_params,
             'time': time_obs,
             'flux': flux_obs,
             'flux_err': flux_err_obs,
             'mags': mag_obs,
             'filt': filts,
             'zp': zp_obs
             }

    return infos

#=========================================================================================================


data = pd.read_parquet("/home/masson/orphans/data/ztf_obsfile_maglimcat.parquet")

data = data.rename(columns={'expMJD': 'mjd', 'filter': 'band', 'fieldID': 'fieldid'})
data['skynoise'] = np.random.normal(size=len(data), loc=50, scale=20)

ztf = skysurvey.ZTF.from_pointings(data)

# take events that happens during the time ZTF was operating
tstart, tstop = ztf.get_timerange()

# load orphan afterglow configurations that are observable by Rubin LSST
orphan_configs_ztf = pd.read_pickle('/home/masson/orphans/data/pseudo_obs/orphan_configs_ztf.pkl')
df_orphans_ztf = pd.DataFrame(data=orphan_configs_ztf)

# simulate observations of orphans by ZTF
obs_ztf_orphans = []

for i in tqdm(range(len(df_orphans_ztf))):

    grb_params = get_grb_params()
    phase, wave, flux = get_flux(grb_params)

    mag_obs = -2.5 * np.log10(np.array(flux) / (2.99792458e18 / (wave ** 2))) - 48.6

    if np.min(mag_obs) < 20.5:
        mag_abs = mag_obs - 5 * np.log10(grb_params['d_L'] / 3.0857e18) + 5

        # =============== #
        #                 #
        #  Orphan         #
        #                 #
        # =============== #

        _ORPHAN_MODEL = get_orphan_model()


        class Orphan(Transient):

            _KIND = "orphan"
            _TEMPLATE = _ORPHAN_MODEL
            _RATE = 10  # event per Gyr**3
            _MODEL = dict(  # when
                t0={"func": np.random.uniform,
                    "kwargs": {"low": tstart, "high": tstop}
                    },

                # what
                magabs={"func": np.random.normal,
                        "kwargs": {"loc": np.min(mag_abs[:, 53]), "scale": 0.}
                        },

                magobs={"func": "magabs_to_magobs",
                        "kwargs": {"z": grb_params['z'], "magabs": "@magabs"}
                        },

                amplitude={"func": "magobs_to_amplitude",
                           "kwargs": {"magobs": "@magobs"}
                           },

                # where
                radec={"func": random_radec,
                       "kwargs": {},
                       "as": ["ra", "dec"]
                       },
            )


        orphan = Orphan()

        oa_transient = Orphan.from_draw(size=10)

        dset = skysurvey.DataSet.from_targets_and_survey(oa_transient, ztf)

        if len(dset.get_ndetection()) > 0 and max(dset.get_ndetection()) >= 5:
            index_max = dset.get_ndetection()[dset.get_ndetection() == max(dset.get_ndetection())].index[0]
            orphan_info = save_ztf_light_curve(grb_params, oa_transient.data.loc[index_max],
                                               dset.data.loc[index_max], index_max)
            obs_ztf_orphans.append(orphan_info.copy())


file = open('/home/masson/orphans/data/pseudo_obs/orphan_pseudo_obs_ztf.pkl', 'wb')
pickle.dump(obs_ztf_orphans, file)
file.close()