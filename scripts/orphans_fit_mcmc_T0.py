import pickle
import emcee
import numpy as np
import afterglowpy as grb

from astropy.cosmology import Planck18 as cosmo
from astropy.time import Time


# ===========================================================================================================================

def flux_to_mag(flux):
    """ Convert flux from milliJansky to AB Magnitude

    1 Jy = 1e-23 erg/cm2/s/Hz
    Fnu = 3631 Jy = 3.631*1e-20 erg/cm2/s/Hz
    ABmag = 0-2.5*log10( Fnu )-48.6 = 0

    :param flux: flux in milli-Jansky
    :return: mag: as the AB Magnitude
    """

    mag = -2.5 * np.log10(flux * 1.0e-26) - 48.6
    return mag


def mag_to_flux(mag):
    """ Convert flux from AB Magnitude to milliJansky

    1 Jy = 1e-23 erg/cm2/s/Hz
    Fnu = 3631 Jy = 3.631*1e-20 erg/cm2/s/Hz
    ABmag = 0-2.5*log10( Fnu )-48.6 = 0

    :param mag: as the AB Magnitude
    :return: flux: flux in milli-Jansky
    """

    flux = pow(10, (26 - (mag + 48.6) / 2.5))
    return flux


def rescale_filters(times, mags, mags_err, filts):
    """
    Calculate the rescaled magnitude to the r-band

    Parameters
    ----------
    times: array
           Time in MJD days
    mags: array
          Observed AB magnitude
    mags_err: array
              Error on the observed magnitude
    filts: array
           Colour of the filter used at the time of the detection: {'u':'b', 'g':'c', 'r':'g', 'i':'orange', 'z':'r', 'Y':'m'}

    Returns
    -------
    time: array
          Times of the rescaled points
    mag_r: array
           Rescaled magnitudes in r-band
    err: array
         Error associated to the observed magnitude
    """

    # colors and mean frequency of the band u, g, r, i, z, y
    filters = ['b', 'c', 'g', 'orange', 'r', 'm']
    all_mean_nu = [840336134453781.4, 629326620516047.8, 482703137570394.2, 397614314115308.1, 344530577088716.56,
                   298760145396604.1]

    filt_obs = filts[filts != 'g']

    if len(filt_obs) != 0:

        unique, counts = np.unique(filt_obs, return_counts=True)
        filt_max = unique[np.argmax(counts)]
        nu_filtmax = all_mean_nu[filters.index(filt_max)]

        mag_r = mags[filts != 'g']
        time_r = times[filts != 'g']

        mag_filtmax = mags[filts == filt_max]
        time_filtmax = times[filts == filt_max]

        flux_r = mag_to_flux(mag_r)
        flux_filtmax = mag_to_flux(mag_filtmax)

        # choose values of -beta between -(p-1)/2 and -p/2
        beta = np.linspace(-0.6, -1.1, 10)

        d = []

        for b in beta:
            # compute the rescaled flux for each beta
            flux_rescaled = flux_filtmax * (all_mean_nu[2] / nu_filtmax) ** (b)

            # compute the euclidean distance between the rescaled flux and the true r-band flux
            d.append(np.sum(np.sqrt((time_filtmax[:, np.newaxis] - time_r[np.newaxis, :]) ** 2 +
                                    (flux_rescaled[:, np.newaxis] - flux_r[np.newaxis, :]) ** 2)))

        beta_min = beta[np.where(d == min(d))]

        all_mag_r = []
        all_err = []
        all_time = []

        # sort the magnitudes, flux, times and errors based on their filter
        for f, nu in zip(filters, all_mean_nu):

            mag_f = mags[filts == f]
            err_f = mags_err[filts == f]
            time_f = times[filts == f]

            flux_f = mag_to_flux(mag_f)

            # for the beta that minimize the distance, rescale the flux for all the bands
            if f != 'g':
                all_mag_r.append(flux_to_mag(flux_f * (all_mean_nu[2] / nu) ** min(beta_min)))
            else:
                all_mag_r.append(flux_to_mag(flux_f))

            all_time.append(time_f)
            all_err.append(err_f)

        # create one array with all the times, one with all the mag in r-band and the associated error
        time = np.concatenate(all_time)
        mag_r = np.concatenate(all_mag_r)
        err = np.concatenate(all_err)

        return time, mag_r, err

    else:

        return times, mags, mags_err


def model(t, nu, params):
    """ Model used to fit, implemented here in the afterglowpy package
    """

    # dictionary containing the fixed parameters
    Z = {'jetType': 4,
         'specType': 0,
         'b': 4,
         'p': 2.2,
         'epsilon_e': 0.1,
         'epsilon_B': 0.01,
         'xi_N': 1.0,
         'z': one_oa['config']['z'],  # redshift is supposed to be known
         'd_L': cosmo.luminosity_distance(one_oa['config']['z']).value * 3.08e24}

    # parameters to fit
    T0, logE, thetaObs, thetaCore, thetaWing, logn = params

    return grb.fluxDensity((t - T0) * grb.day2sec, nu, E0=10 ** logE, thetaObs=thetaObs, thetaCore=thetaCore,
                           thetaWing=thetaWing, n0=10 ** logn, **Z)


def lnlike(p, t, y, yerr):
    """ Compute ln of likelihood
    """

    return -0.5 * np.sum(((y - model(t, nu, p)) / yerr) ** 2)


def lnprior(p):
    """ Compute ln of parameters priors
    """

    T0, logE, thetaObs, thetaCore, thetaWing, logn = p

    # uniform priors
    if (t[0] - 30. < T0 < t[
        0] - 1. and 51. < logE < 55. and 0.0 < thetaCore < np.pi / 2 and thetaCore < thetaWing < np.pi / 2 and 0.0 < thetaObs < np.pi / 2 and np.log10(
            1.) < logn < np.log10(100.)):
        return (0.0)

    return (-np.inf)


def lnpost(p, t, y, yerr):
    """ Compute ln of parameters posteriors
    """

    lp = lnprior(p)

    return lp + lnlike(p, t, y, yerr) if np.isfinite(lp) else -np.inf


def sample_walkers(nsamples, flattened_chain, nu):
    """ Compute the median model and the spread in posteriors of the MCMC
    """

    models = []
    draw = np.floor(np.random.uniform(0, len(flattened_chain), size=nsamples)).astype(int)
    thetas = flattened_chain[draw]

    for i in thetas:
        mod = model(t_fit, nu, i)
        models.append(mod)

    spread = np.std(models, axis=0)
    med_model = np.median(models, axis=0)

    return med_model, spread


# ===========================================================================================================================


all_mean_nu = [840336134453781.4, 629326620516047.8, 482703137570394.2, 397614314115308.1, 344530577088716.56,
               298760145396604.1]

# start = int(sys.argv[1])
# index = [i*100 for i in range(26)]
index = [i * 5 for i in range(520)]
index.append(2602)

file_open = open(f'../data/pseudo_obs/features_po_mvsr_sbat4_5_pts_all_2602.pkl', 'rb')
oa = pickle.load(file_open)
file_open.close()

# sample = oa.iloc[index[start]:index[start + 1]]

oa_dict = oa.to_dict('records')

config = []
grb_time = []
T0_err = []
T0_fit = []
T0_diff = []

for i in range(10,15):
#for i in range(0,10):

    one_oa = oa_dict[i]
    config.append(one_oa['config'])
    grb_time.append(Time(one_oa['grb_time']).mjd)

    print('Rescaling filters...')
    times, mags, mags_err = rescale_filters(np.array(one_oa['time']), np.array(one_oa['mags']),
                                            np.array(one_oa['mags_err']), np.array(one_oa['filt']))
    print('Done!')

    t = times[mags < one_oa['mags_lim']]
    y = mag_to_flux(mags[mags < one_oa['mags_lim']])
    yerr = mag_to_flux(mags[mags < one_oa['mags_lim']] - mags_err[mags < one_oa['mags_lim']]) - y

    print(Time(one_oa['grb_time']).mjd)
    print(t[0])

    # data gathered in a tuple
    data = (t, y, yerr)

    # r-band mean frequency in Hz
    nu = np.array(t.shape)
    nu[:] = all_mean_nu[2]

    z = one_oa['config']['z']

    nwalkers = 30

    if len(t) > 0:
        initial = np.array([t[0] - 10., 52., 0.3, 0.1, 0.2, 0.])
        ndim = len(initial)

        # initial position vector
        p0 = [np.array(initial) + 1e-2 * np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnpost, args=data)
        print(sampler)

        """
        for pos, lnprob, state in sampler.sample(p0, iterations=10):
            best = pos[np.argmax(lnprob), 0]
            print(Time(one_oa['grb_time']).mjd, best)"""

        print("Running burn-in...")
        p0, _, _ = sampler.run_mcmc(p0, 200, progress=True);
        sampler.reset()

        print("Running production...")
        sampler.run_mcmc(p0, 600, progress=True);

        print("Done!")

        samples = sampler.flatchain

        best_params = samples[np.argmax(sampler.flatlnprobability)]
        T0_fit.append(best_params[0])
        print(best_params[0])

T0_err_dict = {'config': np.array(config),
               'T0_fit': np.array(T0_fit),
               'grb_time': np.array(grb_time),
               }

# file = open(f'../data/pseudo_obs/orphans_fit_mcmc_300_test_params_{start}.pkl', 'wb')
file = open(f'../data/pseudo_obs/orphans_fit_mcmc_300.pkl', 'wb')
pickle.dump(T0_err_dict, file)
file.close()
