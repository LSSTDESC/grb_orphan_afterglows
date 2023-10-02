import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pickle
import math
from itertools import chain
import scipy.interpolate


def MinimalMagnitude(lc_open):
    
    """ Save the minimal magnitude (= maximal flux) for each pseudo-observed light curve
    
    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :return: mag_min: list containing the minimal detected magnitude of each configuration
    """
    
    mag_min = []
    
    for lc in lc_open:
        
        # if the afterglow is not in the field of view of Rubin LSST
        if lc == 0:
            mag_min.append(np.nan)

        else:
            # keeping only the observable points, ic values of magnitude that are smaller than the Rubin LSST limiting magnitude
            mag = [lc['mags'][j] for j in range(len(lc['mags'])) if lc['mags'][j] < lc['mags_lim'][j]]
            
            # if no point is observable
            if len(mag) == 0:
                mag_min.append(np.nan)

            # keeping the smallest value of magnitude (brightest point) among all the filters
            else:
                mag_min.append(min(mag))

    return mag_min




def PeakTime(lc_open, data='simu'):

    """ Save the time of the minimal magnitude (= time of the peak) for each pseudo-observed light curve in each filter

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :param data: if the data are from the simulated pseudo-observations ('simu') or from ELAsTiCC ('elasticc') data (default is 'simu')
    :return: time_max: list containing the time of the minimal detected magnitude of each configuration in each filter
    """

    time_max = []
    

    for lc in lc_open:
        
        if lc == 0:
            
            time_max.append(np.nan)
            
        else:
            
            # defining the day of first detection as day zero
            #time_norm = lc_open[i]['time'] - min(lc_open[i]['time'])
            if type(lc['mags_lim']) == float:
                mag = [lc['mags'][j] for j in range(len(lc['mags'])) if lc['mags'][j] < lc['mags_lim']]
            
            else:
                mag = [lc['mags'][j] for j in range(len(lc['mags'])) if lc['mags'][j] < lc['mags_lim'][j]]
                
            
            if len(mag) == 0:
                time_max.append(np.nan)
                
            else:
                time = lc['time']
                #time_norm = lc_open[i]['time'] - min(lc_open[i]['time'])

                if data == 'simu':
                    #time_max.append(time_norm[lc_open[i]['mags'].index(min(lc_open[i]['mags']))])
                    time_max.append(time[lc['mags'].index(min(lc['mags']))])

                elif data == 'elasticc':
                    if math.isnan(lc['mags'][0]) == True:
                        time_max.append(np.nan)
                    else:
                        time_max.append(time[np.where(lc['mags'] == min(lc['mags']))[0]][0])

            # taking the time of the peak in magnitude
            #time_max.append(time_norm[lc_open[i]['mags'].index(min(lc_open[i]['mags']))])
            
    return time_max
    
    
    
def DurationBetweenFirstAndPeak(lc_open):

    """ Save the number of days between the first detection and the peak for each pseudo-observed light curve in each filter

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :return: Dt: list containing the number of days between the first detection and the peak of each configuration in each filter
    """
    
    filterColor = ['b', 'c', 'g', 'orange', 'r', 'm']

    Dt = [[] for _ in range(6)]
    t0 = []
          

    for lc in lc_open:
        
        
        mag_list = []
        time_list = []
        
        if lc == 0:

            Dt = [x.append(np.nan) or x for x in Dt]
                
            
        else:
            
            # sorting each magnitude according to its filter

            for f in filterColor:
    
                mag_list.append([lc['mags'][j] for j in range(len(lc['mags'])) 
                                 if lc['mags'][j] < lc['mags_lim'][j] and lc['filt'][j]==f])
    
                time_list.append([lc['time'][j] for j in range(len(lc['mags'])) 
                                 if lc['mags'][j] < lc['mags_lim'][j] and lc['filt'][j]==f])
        

            all_t0 = []

            for j in range(6):     # for each filter
                
                mag = mag_list[j]
                time = time_list[j]

                if len(mag) != 0:
                    # calculating the difference between the time of first detection and the peak
                    all_t0.append(time[0])
                    Dt[j].append(time[mag.index(min(mag))] - time[0])

                else:
                    Dt[j].append(np.nan)
                    all_t0.append(np.nan)
                
            if np.nanmin(all_t0) != np.nan:
                t0.append(np.nanmin(all_t0))
                
            else:
                t0.append(np.nan)
                    
                    
    return Dt, t0



def Rate(lc_open, data='simu'):

    """ Save the decreasing rates (in the 1/3 and the 3/3 of the decreasing part of the light curve) and the increasing rate for each pseudo-observed light curve in each filter

        :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
        :param data: if the data are from the simulated pseudo-observations ('simu') or from ELAsTiCC ('elasticc') data (default is 'simu')

        Returns
        -------
        rate_dec_1: `list` of `list` of `float`
            a list of the decreasing rate in the 1/3 of the decreasing part of the light curve of each configuration in each filter
        rate_dec_3: `list` of `list` of `float`
            a list of the decreasing rate in the 3/3 of the decreasing part of the light curve of each configuration in each filter
        rate_inc: `list` of `list` of `float`
            a list of the increasing rate of the light curve of each configuration in each filter
        """
    
    
    rate_inc = [[] for _ in range(6)]
    rate_dec_1 = [[] for _ in range(6)]
    rate_dec_3 = [[] for _ in range(6)]
    
    
    for lc in lc_open:
        
        
        if lc == 0:
            
            rate_dec_1 = [x.append(np.nan) or x for x in rate_dec_1]
            rate_dec_3 = [x.append(np.nan) or x for x in rate_dec_1]
            rate_inc = [x.append(np.nan) or x for x in rate_inc]
            
        
        else:
            
            mag_list = []
            time_list = []
            
            
            if data == 'simu':
                
                filterColor = ['b', 'c', 'g', 'orange', 'r', 'm']

                for f in filterColor:

                    mag_list.append([lc['mags'][j] for j in range(len(lc['mags'])) 
                                     if lc['mags'][j] < lc['mags_lim'][j] and lc['filt'][j]==f])

                    time_list.append([lc['time'][j] for j in range(len(lc['mags'])) 
                                      if lc['mags'][j] < lc['mags_lim'][j] and lc['filt'][j]==f])
                    
                    
            elif data == 'elasticc':
                    
                filterColor = ['u', 'g', 'r', 'i', 'z', 'Y']
                
                for f in filterColor:
                                      
                    mag_list.append([lc['mags'][j] for j in range(len(lc['mags'])) 
                                     if lc['mags'][j] < lc['mags_lim'] and lc['filt'][j]==f])

                    time_list.append([lc['time'][j] for j in range(len(lc['mags'])) 
                                      if lc['mags'][j] < lc['mags_lim'] and lc['filt'][j]==f])
                    



            for j in range(6):

                mag = mag_list[j]
                time = time_list[j]

                # if there is at least 4 points and a decreasing phase
                if len(mag) > 3 and mag.index(min(mag)) != len(mag)-1: 

                    mdec = mag[mag.index(min(mag)):len(mag)-1]
                    tdec = time[mag.index(min(mag)):len(mag)-1]


                    if mag.index(min(mag)) != 0:

                        minc = mag[0:mag.index(min(mag))]
                        tinc = time[0:mag.index(min(mag))]
                        dtinc = tinc[minc.index(min(minc))] - tinc[minc.index(max(minc))]
                        
                        if dtinc > 1:
                            rate_inc[j].append((min(minc) - max(minc)) / dtinc)
                            
                        else:
                            rate_inc[j].append(np.nan)

                    else:
                        rate_inc[j].append(np.nan)


                    if len(mdec) > 3:

                        # calculating the rates on the 1/3 and the 3/3 of the light curve
                        m1 = mdec[0:int(len(mdec)/3)]
                        m3 = mdec[2*int(len(mdec)/3):len(mdec)]

                        t1 = tdec[0:int(len(tdec)/3)]
                        t3 = tdec[2*int(len(tdec)/3):len(tdec)]

                        dt1 = t1[m1.index(max(m1))] - t1[m1.index(min(m1))]
                        dt3 = t3[m3.index(max(m3))] - t3[m3.index(min(m3))]

                        # if the 2 points are separated by at least 1 day          
                        if dt1 > 1 and dt3 > 1:
                            rate_dec_1[j].append((max(m1) - min(m1)) / dt1)
                            rate_dec_3[j].append((max(m3) - min(m3)) / dt3)

                        else:                                                                                                                  
                            rate_dec_1[j].append(np.nan)
                            rate_dec_3[j].append(np.nan)                                                                                                             

                    else:

                        rate_dec_1[j].append(np.nan)
                        rate_dec_3[j].append(np.nan)

                elif len(mag) > 3 and mag.index(min(mag)) == len(mag)-1:

                    minc = mag[0:mag.index(min(mag))]
                    tinc = time[0:mag.index(min(mag))]
                    dtinc = tinc[minc.index(min(minc))] - tinc[minc.index(max(minc))]
                    
                    if dtinc > 1:
                        rate_inc[j].append((min(minc) - max(minc)) / dtinc)
                            
                    else:
                        rate_inc[j].append(np.nan)

                    rate_dec_1[j].append(np.nan)
                    rate_dec_3[j].append(np.nan)

                else:
                    rate_inc[j].append(np.nan)
                    rate_dec_1[j].append(np.nan)
                    rate_dec_3[j].append(np.nan)


    return rate_inc, rate_dec_1, rate_dec_3




def Color(lc_open, data='simu'):
    
    """ Calculate the mean g-r color for each pseudo-observed light curve 
    
    expected value for Synchrotron emission = 0.3

    :param lc_open: pseudo-observed light curves save in the lc_configs_*.pkl files
    :param data: if the data are from the simulated pseudo-observations ('simu') or from ELAsTiCC ('elasticc') data (default is 'simu')
    :return: color: list containing the g-r color of each configuration
    """

    color = []    # contains the mean g-r values for each configuration

    for lc in lc_open:
                
        if lc == 0:
            
            color.append(np.nan)
            
        else:
            
            mag_list = []
            time_list = []
            
            
            if data == 'simu':
                        
                filterColor = ['c', 'g']

                for f in filterColor:

                    mag_list.append([lc['mags'][j] for j in range(len(lc['mags'])) 
                                     if lc['mags'][j] < lc['mags_lim'][j] and lc['filt'][j]==f])

                    time_list.append([lc['time'][j] for j in range(len(lc['mags'])) 
                                      if lc['mags'][j] < lc['mags_lim'][j] and lc['filt'][j]==f])
                
            elif data == 'elasticc':
                
                filterColor = ['g', 'r']

                for f in filterColor:

                    mag_list.append([lc['mags'][j] for j in range(len(lc['mags'])) 
                                     if lc['mags'][j] < lc['mags_lim'] and lc['filt'][j]==f])

                    time_list.append([lc['time'][j] for j in range(len(lc['mags'])) 
                                      if lc['mags'][j] < lc['mags_lim'] and lc['filt'][j]==f])

            
            mag_r = mag_list[1]
            mag_g = mag_list[0]


            # removing light curves with only an increasing phase
            if (len(mag_r) != 0 and len(mag_g) != 0 and (mag_r.index(min(mag_r)) != len(mag_r)-1) 
                and (mag_g.index(min(mag_g)) != len(mag_g)-1)):

                time_r = time_list[1]
                time_g = time_list[0]

                # we take points only on the decreasing part of the light curve because the g-r color is constant after the peak
                mdec_r = mag_r[mag_r.index(min(mag_r)):len(mag_r)-1]
                tdec_r = time_r[mag_r.index(min(mag_r)):len(mag_r)-1]

                mdec_g = mag_g[mag_g.index(min(mag_g)):len(mag_g)-1]
                tdec_g = time_g[mag_g.index(min(mag_g)):len(mag_g)-1]


                # we interpolate on the filter with the more points
                # if the r filter has more points, an r point is interpolated
                if (len(mdec_r) > 1 and (tdec_g[len(tdec_g)-1] > tdec_r[0]) and (tdec_r[len(tdec_r)-1] > tdec_g[0])
                    and (tdec_r[len(tdec_r)-1] - tdec_r[0] > tdec_g[len(tdec_g)-1] - tdec_g[0])):

                    interpolation = scipy.interpolate.interp1d(tdec_r, mdec_r)
                    
                    while tdec_g[0] < tdec_r[0]:
                        tdec_g.remove(tdec_g[0])
                        mdec_g.remove(mdec_g[0])
                    while tdec_g[len(tdec_g)-1] > tdec_r[len(tdec_r)-1]:
                        tdec_g.remove(tdec_g[len(tdec_g)-1])
                        mdec_g.remove(mdec_g[len(mdec_g)-1])
                        
                    magr_int = interpolation(tdec_g)
                    color.append(np.mean(mdec_g - magr_int))


                # else, if the g filter has more points, a g point is interpolated
                elif (len(mdec_g) > 1 and (tdec_r[len(tdec_r)-1] > tdec_g[0]) and (tdec_g[len(tdec_g)-1] > tdec_r[0])
                      and (tdec_r[len(tdec_r)-1] - tdec_r[0] < tdec_g[len(tdec_g)-1] - tdec_g[0])):

                    interpolation = scipy.interpolate.interp1d(tdec_g, mdec_g)
                    
                    while tdec_r[0] < tdec_g[0]:
                        tdec_r.remove(tdec_r[0])
                        mdec_r.remove(mdec_r[0])
                    while tdec_r[len(tdec_r)-1] > tdec_g[len(tdec_g)-1]:
                        tdec_r.remove(tdec_r[len(tdec_r)-1])
                        mdec_r.remove(mdec_r[len(mdec_r)-1])
                        
                    magg_int = interpolation(tdec_r)
                    color.append(np.mean(magg_int - mdec_r))


                else:
                    color.append(np.nan)

            else:    
                color.append(np.nan)

                
    return color
    
    
    
def Heatmap(lc_open, mag_min, time_max, t0, rate_dec_1, rate_dec_3, rate_inc, color, Dt, data='simu', parameters='all', annot=True):

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
        
        if data == 'simu':

            # calculating the mean value of (t_max - t0) between each filter
            Dt_one = [Dt[j][i] for j in range(len(Dt)) if Dt[j][i] != np.nan]   
            Dt_mean.append(np.mean(Dt_one))
        
            rate_dec_1_one = [rate_dec_1[j][i] for j in range(len(rate_dec_1)) if math.isnan(rate_dec_1[j][i]) == False] 
            rate_dec_3_one = [rate_dec_3[j][i] for j in range(len(rate_dec_3)) if math.isnan(rate_dec_3[j][i]) == False] 
            rate_inc_one = [rate_inc[j][i] for j in range(len(rate_inc)) if math.isnan(rate_inc[j][i]) == False]

            rate_dec_1_mean.append(np.mean(rate_dec_1_one))
            rate_dec_3_mean.append(np.mean(rate_dec_3_one))
            rate_inc_mean.append(np.mean(rate_inc_one))
            
        elif data == 'elasticc':
            
            rate_dec_1_one = [rate_dec_1[j][i] for j in range(len(rate_dec_1)) if math.isnan(rate_dec_1[j][i]) == False] 
            rate_dec_3_one = [rate_dec_3[j][i] for j in range(len(rate_dec_3)) if math.isnan(rate_dec_3[j][i]) == False] 
            rate_inc_one = [rate_inc[j][i] for j in range(len(rate_inc)) if math.isnan(rate_inc[j][i]) == False]

            rate_dec_1_mean.append(np.mean(rate_dec_1_one))
            rate_dec_3_mean.append(np.mean(rate_dec_3_one))
            rate_inc_mean.append(np.mean(rate_inc_one))
            

    
    results = pd.DataFrame(list(zip(mag_min, time_max, t0, rate_dec_1_mean, rate_dec_3_mean, rate_inc_mean, color, Dt_mean)), columns=['mag_peak', 'time_peak', 'first_detect', 'rate_dec_1', 'rate_dec_3', 'rate_inc', 'color', 'dt'])
    
    results_not_zeros = results[(results['time_peak'] > 0) & (results['dt'] > 0)] 
    
    
    if data == 'simu':
    
        lc_param = []

        dict_nan = {'jetType': -1,
              'specType': 0,
              'b': 4,
              'thetaObs': np.nan,
              'E0': np.nan,
              'thetaWing': np.nan,
              'thetaCore': np.nan,
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

        df = pd.DataFrame(lc_param, columns = ['E0', 'n0', 'z', 'thetaObs', 'thetaCore', 'thetaWing'])

        sns.set_theme(style="whitegrid")

        df['log10(E0)'] = np.log10(df['E0'])
        df['log10(n0)'] = np.log10(df['n0'])
        df['cos(thetaObs)'] = np.cos(df['thetaObs'])

        all_results = pd.concat([df, results], axis=1)
        display(all_results.dropna(axis=0))
        #all_results['log(time_max)'] = np.log10(all_results['time_max'])
        #all_results['log(tmax-t0)'] = np.log10(all_results['tmax-t0'])

        all_results_not_zeros = pd.concat([df, results_not_zeros], axis=1)
        #all_results_not_zeros['log(time_max)'] = np.log10(all_results_not_zeros['time_max'])
        #all_results_not_zeros['log(tmax-t0)'] = np.log10(all_results_not_zeros['tmax-t0'])
        
        
    elif data == 'elasticc':
        
        all_results_not_zeros = results

    
    plt.figure(figsize=(15, 13))
    
    
    if parameters == 'model':
        correlations = df.corr(numeric_only=True)
        sns.set(font_scale=1.7)
        sns.heatmap(df.corr(numeric_only=True), annot=annot, annot_kws={"size": 18}, center=0)
        plt.xticks(rotation=90) 
        plt.yticks(rotation=0) 
        g = sns.PairGrid(df.corr(numeric_only=True), diag_sharey=False, corner=True)
        g.map_upper(sns.histplot)
        g.map_lower(sns.kdeplot, fill=True)
        g.map_diag(sns.kdeplot, color='midnightblue', lw=1.5)
        g.map_diag(sns.histplot, kde=True);
        
    elif parameters == 'all':
        correlations = all_results_not_zeros.corr(numeric_only=True)
        sns.set(font_scale=1.3)
        sns.heatmap(all_results_not_zeros.corr(numeric_only=True), annot=annot, annot_kws={"size": 11}, center=0)
        plt.xticks(rotation=90) 
        plt.yticks(rotation=0) 
        g = sns.PairGrid(all_results_not_zeros.corr(numeric_only=True), diag_sharey=False, corner=True)
        g.map_lower(sns.kdeplot, fill=True)
        g.map_diag(sns.kdeplot, lw=1.5)
        #g.map_diag(sns.histplot, kde=True);
    
    elif parameters == 'features':
        correlations = results_not_zeros.corr(numeric_only=True)
        sns.set(font_scale=1.7)
        sns.heatmap(results_not_zeros.corr(numeric_only=True), annot=annot, annot_kws={"size": 18}, center=0)
        plt.xticks(rotation=90) 
        plt.yticks(rotation=0) 
        g = sns.PairGrid(results_not_zeros.corr(numeric_only=True), diag_sharey=False, corner=True)
        g.map_lower(sns.kdeplot, fill=True)
        g.map_diag(sns.kdeplot, lw=1.5)
        #g.map_diag(sns.histplot, kde=True);
    
    return correlations
