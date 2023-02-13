import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pickle
import math
from itertools import chain


def MinimalMagnitude(lc_open):
    
    """ Save the minimal magnitude (= maximal flux) for each pseudo-observed light curve
    
    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :return: mag_min: list containing the minimal detected magnitude of each configuration
    """
    
    mag_min = []
    
    for i in range(len(lc_open)):
            
        if lc_open[i] == 0:     # if the afterglow is not observed by Rubin LSST
            mag_min.append(np.nan)
		    
        else:
            # keeping only the values of magnitude that are smaller than the Rubin LSST limiting magnitude
            mag = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j]]
            
		# if no point is observable
        if len(mag) == 0:
            mag_min.append(np.nan)
            
        # keeping the smallest value of magnitude among all the filters
        else:
            mag_min.append(min(mag))
		
    return mag_min
	
	
	
	
def PeakTime(lc_open):

    """ Save the time of the minimal magnitude (= time of the peak) for each pseudo-observed light curve in each filter

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :return: time_max: list containing the time of the minimal detected magnitude of each configuration in each filter
    """

    time_max = []

    for i in range(len(lc_open)):
        
        if lc_open[i] == 0:
            
            time_max.append(np.nan)
            
        else:
            
            # defining the day of first detection as day zero
            time_norm = lc_open[i]['time'] - min(lc_open[i]['time'])
            
            # taking the time of the peak in magnitude
            time_max.append(time_norm[lc_open[i]['mags'].index(min(lc_open[i]['mags']))])
            
    return time_max
    
    
    
def DurationBetweenFirstAndPeak(lc_open):

    """ Save the number of days between the first detection and the peak for each pseudo-observed light curve in each filter

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :return: Dt: list containing the number of days between the first detection and the peak of each configuration in each filter
    """

    Dt = [[] for _ in range(6)]


    for i in range(len(lc_open)):
        
        if lc_open[i] == 0:
            
            for j in range(6):

                Dt[j].append(np.nan)
                
            
        else:
            
            # sorting each magnitude according to its filter

            mag_u = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='b']
            time_u = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='b']

            mag_g = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='c']
            time_g = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='c']

            mag_r = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='g']
            time_r = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='g']

            mag_i = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='orange']
            time_i = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='orange']

            mag_z = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='r']
            time_z = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='r']

            mag_y = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='m']
            time_y = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='m']

            mag_list = [mag_u, mag_g, mag_r, mag_i, mag_z, mag_y]
            time_list = [time_u, time_g, time_r, time_i, time_z, time_y]

                
            for j in range(6):     # for each filter
                
                mag = mag_list[j]
                
                if len(mag) != 0:
                    
                    # calculating the difference between the time of first detection and the peak
                    Dt[j].append(time_list[j][mag.index(min(mag))] - time_list[j][0])
                
                else:
                    Dt[j].append(np.nan)
                    
    return Dt




def Rate(lc_open):

    """ Save the decreasing rates (in the 1/3 and the 3/3 of the decreasing part of the light curve) and the increasing rate for each pseudo-observed light curve in each filter

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    
    Returns
    -------
    rate_dec_1: `list` of `list` of `float`
        a list of the decreasing rate in the 1/3 of the decreasing part of the light curve of each configuration in each filter
    rate_dec_3: `list` of `list` of `float`
        a list of the decreasing rate in the 3/3 of the decreasing part of the light curve of each configuration in each filter
    rate_inc: `list` of `list` of `float`
        a list of the increasing rate of the light curve of each configuration in each filter
    """


    rate_dec_1 = [[] for _ in range(6)]    # decreasing rate of the flux in the first third of the curve (< 0 because the magnitude increases)
    rate_dec_3 = [[] for _ in range(6)]    # decreasing rate of the flux in the third third of the curve

    rate_inc = [[] for _ in range(6)]    # increasing rate of the flux (> 0 because the magnitude decreases)


    for i in range(len(lc_open)):
        
        if lc_open[i] == 0:
            
            for j in range(6):

                rate_dec_1[j].append(np.nan)
                rate_dec_3[j].append(np.nan)
                
                rate_inc[j].append(np.nan)
            
        else:

            mag_u = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='b']
            time_u = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='b']

            mag_g = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='c']
            time_g = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='c']

            mag_r = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='g']
            time_r = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='g']

            mag_i = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='orange']
            time_i = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='orange']

            mag_z = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='r']
            time_z = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='r']

            mag_y = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='m']
            time_y = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='m']

            mag_list = [mag_u, mag_g, mag_r, mag_i, mag_z, mag_y]
            time_list = [time_u, time_g, time_r, time_i, time_z, time_y]

                
            for j in range(6):
                
                mag = mag_list[j]
                time = time_list[j]
                
                # if there is at least one point and if the minimal magnitude happens before the maximal magnitude
                if len(mag) > 3 and mag.index(min(mag)) < mag.index(max(mag)):        

                    # if we only have the decreasing phase
                    if mag.index(min(mag)) == 0:
                        
                        mag_1 = mag[0:int(len(mag)/3)]
                        mag_3 = mag[2*int(len(mag)/3):len(mag)]

                        # calculation of the decreasing rate as : rate_dec = (mag_max - min_max) / dt 
                        if (time[mag_1.index(max(mag_1))] - time[mag_1.index(min(mag_1))]) > 1 and (time[mag_3.index(max(mag_3))] - time[mag_3.index(min(mag_3))]) > 1:
                            
                            rate_dec_1[j].append((mag_1[mag_1.index(max(mag_1))] - mag_1[mag_1.index(min(mag_1))]) / (time[mag_1.index(max(mag_1))] - time[mag_1.index(min(mag_1))]))
                            rate_dec_3[j].append((mag_3[mag_3.index(max(mag_3))] - mag_3[mag_3.index(min(mag_3))]) / (time[mag_3.index(max(mag_3))] - time[mag_3.index(min(mag_3))]))
                            
                        else:
                                                                                                                                           
                            rate_dec_1[j].append(np.nan)
                            rate_dec_3[j].append(np.nan)
                                                                                                                                       
                        rate_inc[j].append(np.nan)
            
                    # else, if we have a increasing and decreasing phase
                    else:

                        mag_inc = mag[0:mag.index(min(mag))]
                        mag_dec = mag[mag.index(min(mag)):len(mag)-1]
                        
                        rate_inc[j].append((mag[mag.index(min(mag))] - mag[mag.index(max(mag_inc))]) / (time[mag.index(min(mag))] - time[mag.index(max(mag_inc))]))
                        
                        if len(mag_dec) > 3:
                        
                            mag_1 = mag_dec[0:int(len(mag_dec)/3)]
                            mag_3 = mag_dec[2*int(len(mag_dec)/3):len(mag_dec)]
                            
                            if (time[mag_1.index(max(mag_1))] - time[mag_1.index(min(mag_1))]) > 1 and (time[mag_3.index(max(mag_3))] - time[mag_3.index(min(mag_3))]) > 1:

                                rate_dec_1[j].append((mag_1[mag_1.index(max(mag_1))] - mag_1[mag_1.index(min(mag_1))]) / (time[mag_1.index(max(mag_1))] - time[mag_1.index(min(mag_1))]))                                
                                rate_dec_3[j].append((mag_3[mag_3.index(max(mag_3))] - mag_3[mag_3.index(min(mag_3))]) / (time[mag_3.index(max(mag_3))] - time[mag_3.index(min(mag_3))]))
     
                            else:
                                                                                                                                           
                                rate_dec_1[j].append(np.nan)
                                rate_dec_3[j].append(np.nan)                                                                                                           
                                                                                                                                           
                        else:
                            
                            rate_dec_1[j].append(np.nan)
                            rate_dec_3[j].append(np.nan)
                        
                    
                    
                # if there is at least one point and if the minimal magnitude happens after the maximal magnitude
                elif len(mag) > 3 and mag.index(min(mag)) > mag.index(max(mag)):

                    # if we only have the increasing phase
                    if mag.index(min(mag)) == (len(mag)-1):

                        rate_inc[j].append((mag[mag.index(min(mag))] - mag[mag.index(max(mag))]) / (time[mag.index(min(mag))] - time[mag.index(max(mag))]))

                        rate_dec_1[j].append(np.nan)
                        rate_dec_3[j].append(np.nan)

                    # else, if we have a increasing and decreasing phase
                    else:
                        
                        mag_inc = mag[0:mag.index(min(mag))]
                        mag_dec = mag[mag.index(min(mag)):len(mag)-1]

                        rate_inc[j].append((mag[mag.index(min(mag))] - mag[mag.index(max(mag))]) / (time[mag.index(min(mag))] - time[mag.index(max(mag))]))

                        if len(mag_dec) > 3:
                        
                            mag_1 = mag_dec[0:int(len(mag_dec)/3)]
                            mag_3 = mag_dec[2*int(len(mag_dec)/3):len(mag_dec)]
                            
                        
                            if (time[mag_1.index(max(mag_1))] - time[mag_1.index(min(mag_1))]) > 1 and (time[mag_3.index(max(mag_3))] - time[mag_3.index(min(mag_3))]) > 1:

                                rate_dec_1[j].append((mag_1[mag_1.index(max(mag_1))] - mag_1[mag_1.index(min(mag_1))]) / (time[mag_1.index(max(mag_1))] - time[mag_1.index(min(mag_1))]))
                                rate_dec_3[j].append((mag_3[mag_3.index(max(mag_3))] - mag_3[mag_3.index(min(mag_3))]) / (time[mag_3.index(max(mag_3))] - time[mag_3.index(min(mag_3))]))
                        
                            else:
                                                                                                                                           
                                rate_dec_1[j].append(np.nan)
                                rate_dec_3[j].append(np.nan)                                                                                                               
                                                                                                                                           
                        else:
                            
                            rate_dec_1[j].append(np.nan)
                            rate_dec_3[j].append(np.nan)
                
                else:
                    
                    rate_inc[j].append(np.nan)
                    
                    rate_dec_1[j].append(np.nan)
                    rate_dec_3[j].append(np.nan)
                    
    return rate_dec_1, rate_dec_3, rate_inc
    
    


def Color(lc_open):

    """ Calculate the mean g-r color for each pseudo-observed light curve 
    
    expected value for Synchrotron emission = 0.3

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :return: color: list containing the g-r color of each configuration
    """

    color = []    # contains the mean g-r values for each configuration

    for i in range(len(lc_open)):
        
        color_one = []     # contains all the g-r values for one configuration
        
        if lc_open[i] == 0:
            
            color.append(np.nan)
            
        else:
        
            mag_g = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='c']
            time_g = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='c']

            mag_r = [lc_open[i]['mags'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='g']
            time_r = [lc_open[i]['time'][j] for j in range(len(lc_open[i]['mags'])) if lc_open[i]['mags'][j] < lc_open[i]['mags_lim'][j] and lc_open[i]['filt'][j]=='g']

        
            if len(mag_r) != 0 and len(mag_g) != 0:

                # we take points only after the max because the g-r color is constant after the peak
                
                if time_r[0] <= time_g[0]:
                    
                    time_r_pm = [time_r[i] for i in range(len(time_r)) if time_r[i] > time_r[mag_r.index(min(mag_r))]]
                    mag_r_pm = [mag_r[time_r.index(time_r_pm[i])] for i in range(len(time_r_pm))]
                    time_g_pm = [time_g[i] for i in range(len(time_g)) if time_g[i] > time_r[mag_r.index(min(mag_r))]]
                    mag_g_pm = [mag_g[time_g.index(time_g_pm[i])] for i in range(len(time_g_pm))]
                    
                else:
                    
                    time_g_pm = [time_g[i] for i in range(len(time_g)) if time_g[i] > time_g[mag_g.index(min(mag_g))]]
                    mag_g_pm = [mag_g[time_g.index(time_g_pm[i])] for i in range(len(time_g_pm))]
                    time_r_pm = [time_r[i] for i in range(len(time_r)) if time_r[i] > time_g[mag_g.index(min(mag_g))]]
                    mag_r_pm = [mag_r[time_r.index(time_r_pm[i])] for i in range(len(time_r_pm))]
                
                
                # we extrapolate on the filter with the more points
                # if the r filter has more points
                if len(mag_r_pm) != 0 and len(mag_g_pm) != 0 and len(mag_r_pm) > len(mag_g_pm):

                    for j in range(len(mag_r_pm)):

                        for k in range(len(mag_g_pm)):
                            
                            # if the g and r points are less than half-day apart, if the g point is after the r point and not the last point
                            if (time_g_pm[k] - time_r_pm[j] < 0.5) and (time_g_pm[k] - time_r_pm[j] > 0) and j != len(time_r_pm)-1:

                                # linear extrapolation of the r point if it was taken at the same moment as the g point
                                # the g point is taken after the r point so the extrapolation is made between the r point i and i+1
                                a = (mag_r_pm[j+1] - mag_r_pm[j]) / (time_r_pm[j+1] - time_r_pm[j])
                                mag_r_int = a * (time_g_pm[k] - time_r_pm[j]) + mag_r_pm[j]

                                color_one.append(mag_g_pm[k] - mag_r_int)

                            # same but the g point is before the r point
                            elif (time_g_pm[k] - time_r_pm[j] > -0.5) and (time_g_pm[k] - time_r_pm[j] < 0) and j != 0:

                                a = (mag_r_pm[j-1] - mag_r_pm[j]) / (time_r_pm[j-1] - time_r_pm[j])
                                mag_r_int = a * (time_g_pm[k] - time_r_pm[j]) + mag_r_pm[j]

                                color_one.append(mag_g_pm[k] - mag_r_int)
                  
                
                # if the g filter has more points                
                elif len(mag_r_pm) != 0 and len(mag_g_pm) != 0 and len(mag_r_pm) < len(mag_g_pm):
                    
                    for j in range(len(mag_g_pm)):

                        for k in range(len(mag_r_pm)):

                            if (time_r_pm[k] - time_g_pm[j] < 0.5) and (time_r_pm[k] - time_g_pm[j] > 0) and j != len(time_g_pm)-1:

                                a = (mag_g_pm[j+1] - mag_g_pm[j]) / (time_g_pm[j+1] - time_g_pm[j])
                                mag_g_int = a * (time_r_pm[k] - time_g_pm[j]) + mag_g_pm[j]

                                color_one.append(mag_g_int - mag_r_pm[k])

                            elif (time_r_pm[k] - time_g_pm[j] > -0.5) and (time_r_pm[k] - time_g_pm[j] < 0) and j != 0:

                                a = (mag_g_pm[j-1] - mag_g_pm[j]) / (time_g_pm[j-1] - time_g_pm[j])
                                mag_g_int = a * (time_r_pm[k] - time_g_pm[j]) + mag_g_pm[j]

                                color_one.append(mag_g_int - mag_r_pm[k])
                
                # we save the mean of all the g-r values for one light curve to have one g-r value for each configuration
                color.append(np.mean(color_one))
                
                
            else:
                color.append(np.nan)
                
    return color
    
    
    
    
def Heatmap(lc_open, mag_min, time_max, rate_dec_1, rate_dec_3, rate_inc, color, Dt, parameters='all', annot=True):

    """ Calculate the correlations and plot the correlation matrix between the model parameters and the features

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :param mag_min: `list` from MinimalMagnitude function
    :param time_max: `list` from PeakTime function
    :param rate_dec_1, rate_dec_3, rate_inc: `list` from Rates function
    :param color: `list` from Color function
    :param Dt: `list` from DurationBetweenFirstAndPeak function
    
    :return: correlations : panda dataframe containing the correlation between each parameter
    """

    Dt_mean = []
    
    rate_dec_1_mean = []
    rate_dec_3_mean = []
    rate_inc_mean = []


    for i in range(len(lc_open)):

        # calculating the mean value of (t_max - t0) between each filter
        Dt_one = [Dt[j][i] for j in range(len(Dt)) if Dt[j][i] != np.nan]   
        Dt_mean.append(np.mean(Dt_one))
        
        rate_dec_1_one = [rate_dec_1[j][i] for j in range(len(rate_dec_1)) if math.isnan(rate_dec_1[j][i]) == False] 
        rate_dec_3_one = [rate_dec_3[j][i] for j in range(len(rate_dec_3)) if math.isnan(rate_dec_3[j][i]) == False] 
        rate_inc_one = [rate_inc[j][i] for j in range(len(rate_inc)) if math.isnan(rate_inc[j][i]) == False]
        
        rate_dec_1_mean.append(np.mean(rate_dec_1_one))
        rate_dec_3_mean.append(np.mean(rate_dec_3_one))
        rate_inc_mean.append(np.mean(rate_inc_one))

    
    results = pd.DataFrame(list(zip(mag_min, time_max, rate_dec_1_mean, rate_dec_3_mean, rate_inc_mean, color, Dt_mean)), columns=['mag_min', 'time_max', 'rate_dec_1', 'rate_dec_3', 'rate_inc', 'color', 'tmax-t0'])
    
    results_not_zeros = results[(results['time_max'] > 0) & (results['tmax-t0'] > 0)] 
    
    lc_param = []

    dict_nan = {'jetType': -1,
          'specType': 0,
          'b': 4,
          'thetaObs': np.nan,
          'E0': np.nan,
          'thetaWing': np.nan,
          'thetaCore': 0.15,
          'n0': np.nan,
          'p': 2.2,
          'epsilon_e': 0.1,
          'epsilon_B': 0.01,
          'xi_N': 1.0,
          'd_L': 1.0e+28,
          'z': np.nan}


    for i in range(len(lc_open)):
        if lc_open[i] == 0:
            lc_param.append(dict_nan)
        else:
            lc_param.append(lc_open[i]['config'])
            
    df = pd.DataFrame(lc_param, columns = ['E0', 'n0', 'z', 'thetaObs', 'thetaWing'])
    
    sns.set_theme(style="whitegrid")

    df['log10(E0)'] = np.log10(df['E0'])
    df['log10(n0)'] = np.log10(df['n0'])
    df['cos(thetaObs)'] = np.cos(df['thetaObs'])
    
    all_results = pd.concat([df, results], axis=1)
    all_results['log(time_max)'] = np.log10(all_results['time_max'])
    all_results['log(tmax-t0)'] = np.log10(all_results['tmax-t0'])
    
    all_results_not_zeros = pd.concat([df, results_not_zeros], axis=1)
    all_results_not_zeros['log(time_max)'] = np.log10(all_results_not_zeros['time_max'])
    all_results_not_zeros['log(tmax-t0)'] = np.log10(all_results_not_zeros['tmax-t0'])
    
    plt.figure(figsize=(15, 13))
    
    if parameters == 'model':
        correlations = df.corr()
        sns.heatmap(df.corr(), annot=annot, center=0);
        
    elif parameters == 'all':
        correlations = all_results.corr()
        sns.heatmap(all_results.corr(), annot=annot, center=0);
    
    return correlations
