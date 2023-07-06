"""Plotting afterglow light curve module

"""

import numpy as np
import afterglowpy as grb
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from astropy.time import Time


from orphans.tools import ObsTime, mag_to_flux
from orphans.tools import ObsTime, mag_to_flux, galactic_extinction



def plot_simulation(t, fnu, flux='mag'):

    """ Plot theoretical light curve of an orphan afterglow

    :param t: time for which the flux is calculated
    :param fnu: value of the flux of the afterglow in mJy
    :param flux: whether you want to plot the flux in magnitude or in mJy. Takes the values 'mag' or 'flux', default is 'mag'

    :return: plot of the afterglow light curve
    """

    plt.figure(figsize=(15, 12))
    plt.rcParams.update({'font.size': 22})

    if flux == 'mag':
        mag = -2.5 * np.log10(fnu * 1.0e-26) - 48.6
        plt.plot(t, mag, color='royalblue')
        plt.ylabel('AB Magnitude')
        plt.gca().invert_yaxis()
        plt.ylim(30, 15)
        plt.axhline(y=24.5, color='black', linewidth=0.5, label='Rubin/LSST nightly limiting magnitude')

    elif flux == 'flux':
        plt.plot(t, fnu, color='royalblue')
        plt.ylabel('Fux (mJy)')
        plt.yscale('log')
        plt.ylim(mag_to_flux(30), mag_to_flux(15))
        plt.axhline(y=mag_to_flux(24.5), color='black', linewidth=0.5,
                    label='Rubin/LSST nightly limiting magnitude')

    plt.xlabel('Time (days)')
    plt.xscale('log')

    plt.show()




def plot_pseudo_obs(lc, path_dustmaps, lc_th=True, flux='mag', extinction=True):

    """ Plot pseudo-observed light curve of an afterglow

    :param lc: dictionary containing the pseudo-observed light curve, calculated with the generate_pseudo_obs function (in the pickling module)
    :param path_dustmaps: path of the directory where you the conversion table is saved
    :param flux: whether you want to plot the flux in magnitude or in mJy. Takes the values 'mag' or 'flux', default is 'mag'
    :param extinction: whether you want to consider the galactic extinction. Takes the values 'True' or 'False', default is 'True'

    :return: plot of the afterglow pseudo-observed light curve
    """

    #plt.rcParams.update({'font.size': 18})
    #plt.rcParams["figure.figsize"] = [15, 12]
    
    plt.rcParams.update({'font.size': 35})
    plt.figure(figsize=(20, 15))

    filterlist = ['u', 'g', 'r', 'i', 'z', 'y']
    filtercolors = {'u': 'b', 'g': 'c', 'r': 'g', 'i': 'orange', 'z': 'r', 'y': 'm'}
    all_nu = [840336134453781.4, 629326620516047.8, 482703137570394.2, 397614314115308.1, 344530577088716.56, 298760145396604.1]
    
    x_times = lc['time']
    z_colors = lc['filt']

    if flux == 'mag':
        y_mags = lc['mags']
        mags_lim = lc['mags_lim']
        mags_err = lc['mags_err']
        plt.gca().invert_yaxis()
        plt.ylabel('Observed magnitude', fontsize=35)
    elif flux == 'flux':
        y_mags = list(map(mag_to_flux, lc['mags']))
        mags_lim = list(map(mag_to_flux, lc['mags_lim']))
        mags_err = list(map(mag_to_flux, lc['mags_err']))
        plt.ylabel('Observed flux (mJy)', fontsize=35)
        plt.yscale('log')

    if lc_th == True:
        
        t = np.geomspace(min(x_times)-Time(lc['grb_time']).mjd, max(x_times)-Time(lc['grb_time']).mjd, num=100)

        for nu in all_nu:
            
            if extinction == True:
                a_lambda_u, a_lambda_g, a_lambda_r, a_lambda_i, a_lambda_z, a_lambda_y = galactic_extinction(lc['grb_coord'], path_dustmaps)
                a_lambda = [a_lambda_u, a_lambda_g, a_lambda_r, a_lambda_i, a_lambda_z, a_lambda_y]
            
            else: 
                a_lambda = [0, 0, 0, 0, 0, 0]
            
            mag = -2.5 * np.log10(grb.fluxDensity(t*grb.day2sec, nu, **lc['config'])*1.0e-26) - 48.6 + a_lambda[all_nu.index(nu)]
            if nu == all_nu[0]:
                plt.plot(t+Time(lc['grb_time']).mjd, mag, color=filtercolors[filterlist[all_nu.index(nu)]], alpha=0.3, label='model')
            else:
                plt.plot(t+Time(lc['grb_time']).mjd, mag, color=filtercolors[filterlist[all_nu.index(nu)]], alpha=0.3)

    # plot pseudo observed light curve
    for x, y, z, m, e in zip(x_times, y_mags, z_colors, mags_lim, mags_err):
        if flux == 'mag' and y < m:
            plt.scatter(x, y, c=z, s=100)
            plt.errorbar(x, y, e, c=z, capsize=4)
            plt.scatter(x, m, c=z, marker='2', s=100)
        elif flux == 'flux' and y > m:
            plt.scatter(x, y, c=z, s=100)
            plt.errorbar(x, y, e, c=z, capsize=4)
            plt.scatter(x, m, c=z, marker='2', s=100)

    #plt.title('Pseudo observed light curve', fontsize=28)
    plt.xlabel('Time in MJD', fontsize=35)

    legend_elements_2 = list()

    for filt in filterlist:
        fcolor = filtercolors[filt]
        legend_elements_2.append(Line2D([0], [0], marker='o', color=fcolor, label=filt,
                                        markerfacecolor=fcolor, markersize=10))
    plt.legend(handles=legend_elements_2)
    #plt.legend()
    
    #image_format = 'png' # e.g .png, .svg, etc.
    #image_name = 'pseudo_ob_example.png'

    #plt.savefig(image_name, format=image_format, dpi=1200)
    plt.show()