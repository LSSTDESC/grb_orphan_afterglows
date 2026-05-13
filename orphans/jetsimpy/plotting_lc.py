"""Plot afterglow light curves module
"""

import numpy as np
import jetsimpy
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from astropy.time import Time

from modules.tools import flux_to_mag


def plot_pseudo_obs(lc, lc_th=True):

    """ Plot pseudo-observed light curve of an orphan afterglow

    :param lc: dictionary containing the pseudo-observed light curve, calculated with the generate_pseudo_obs function (in the pickling module)

    :return: plot of the afterglow pseudo-observed light curve
    """


    plt.rcParams.update({'font.size': 18})
    plt.figure(figsize=(10,8))

    filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
    filtercolors = {'u': 'b', 'g': 'c', 'r': 'g', 'i': 'orange', 'z': 'r', 'y': 'm'}
    all_nu = [840336134453781.4, 629326620516047.8, 482703137570394.2, 397614314115308.1, 344530577088716.56, 298760145396604.1]
    
    x_times = lc['time']
    z_colors = lc['filt']
    y_mags = lc['mags']
    mags_lim = lc['mags_lim']
    mags_err = lc['mags_err']

    P = lc['config']


    if lc_th == True:
        
        t_min = np.min(np.array(x_times)[np.array(y_mags) < np.array(mags_lim)])
        t_max = np.max(np.array(x_times)[np.array(y_mags) < np.array(mags_lim)])
        t = np.geomspace(t_min-Time(lc['grb_time']).mjd, t_max-Time(lc['grb_time']).mjd, num=100)

        tsecond = t*3600*24

        for nu in all_nu:

            jet = jetsimpy.Jet(
                jetsimpy.PowerLaw(P["theta_c"], P["Eiso"], lf0=P["lf"]),  # jet profile
                P["A"],  # wind number density scale
                P["n0"],  # ism number density scale
                spread=False,  # w/wo spreading effect
                grid=jetsimpy.ForwardJetRes(P["theta_c"], 129)  # resolution
            )

            flux = jet.FluxDensity(
                tsecond,  # [second] observing time span
                nu,  # [Hz]     observing frequency
                P,  # parameter dictionary
            )

            mag = flux_to_mag(flux)

            if nu == all_nu[0]:
                plt.plot(t+Time(lc['grb_time']).mjd, mag, color=filtercolors[filterlist[all_nu.index(nu)]], alpha=0.3, label='model')
            else:
                plt.plot(t+Time(lc['grb_time']).mjd, mag, color=filtercolors[filterlist[all_nu.index(nu)]], alpha=0.3)

    # plot pseudo observed light curve
    for x, y, z, m, e in zip(x_times, y_mags, z_colors, mags_lim, mags_err):
        if y < m:
            plt.scatter(x, y, c=z, s=20)
            plt.errorbar(x, y, e, c=z, capsize=0)
            plt.scatter(x, m, c=z, marker='2', s=60)


    #plt.title('Pseudo observed light curve', fontsize=28)
    plt.xlabel('Time (MJD)', fontsize=20)
    plt.gca().invert_yaxis()
    plt.ylabel('AB Magnitude', fontsize=20)

    legend_elements_2 = list()

    for filt in filterlist:
        fcolor = filtercolors[filt]
        legend_elements_2.append(Line2D([0], [0], marker='o', color=fcolor, label=filt,
                                        markerfacecolor=fcolor, markersize=5))
    plt.legend(handles=legend_elements_2)

    plt.show()