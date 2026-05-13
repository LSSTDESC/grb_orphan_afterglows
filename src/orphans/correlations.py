"""
Utility functions for computing observational features from pseudo-observed light curves.
These functions are used by notebooks and the documentation. Public API functions are
provided in snake_case, with deprecated CamelCase aliases for backward compatibility.
"""

import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.interpolate
from IPython.display import display

# ---------------------------------------------------------------------------
# Helper aliases (deprecated) – simple name binding, no wrapper function.
# ---------------------------------------------------------------------------

# The aliases are assigned after the real functions are defined.

# ---------------------------------------------------------------------------
# Core functions – snake_case API
# ---------------------------------------------------------------------------

def minimal_magnitude(lc_open):
    """Return the minimal magnitude (maximal flux) for each pseudo-observed light curve.

    Parameters
    ----------
    lc_open : iterable
        Pseudo-observed light curves saved in the ``lc_configs_*.pkl`` files.

    Returns
    -------
    list
        Minimal detected magnitude for each configuration (``np.nan`` when not
        observable).
    """

    mag_min = []
    for lc in lc_open:
        if lc == 0:
            mag_min.append(np.nan)
        else:
            # Keep only observable points (magnitudes brighter than the limit).
            mag = [
                lc["mags"][j]
                for j in range(len(lc["mags"]))
                if lc["mags"][j] < lc["mags_lim"][j]
            ]
            mag_min.append(min(mag) if mag else np.nan)
    return mag_min


def peak_time(lc_open, data="simu"):
    """Return the time of the minimal magnitude (peak) for each pseudo-observed light curve.

    Parameters
    ----------
    lc_open : iterable
        Pseudo-observed light curves.
    data : str, optional
        ``"simu"`` for simulated pseudo-observations (default) or ``"elasticc"``.

    Returns
    -------
    list
        Time of the minimal detected magnitude for each configuration (``np.nan``
        when not available).
    """

    time_max = []
    for lc in lc_open:
        if lc == 0:
            time_max.append(np.nan)
        else:
            if isinstance(lc["mags_lim"], float):
                mag = [
                    lc["mags"][j]
                    for j in range(len(lc["mags"]))
                    if lc["mags"][j] < lc["mags_lim"]
                ]
            else:
                mag = [
                    lc["mags"][j]
                    for j in range(len(lc["mags"]))
                    if lc["mags"][j] < lc["mags_lim"][j]
                ]
            if not mag:
                time_max.append(np.nan)
            else:
                time = lc["time"]
                if data == "simu":
                    time_max.append(time[lc["mags"].index(min(lc["mags"]))])
                elif data == "elasticc":
                    if math.isnan(lc["mags"][0]):
                        time_max.append(np.nan)
                    else:
                        time_max.append(
                            time[np.where(lc["mags"] == min(lc["mags"]))[0]][0]
                        )
    return time_max


def duration_between_first_and_peak(lc_open):
    """Return days between first detection and peak for each pseudo-observed
    light curve.

    Parameters
    ----------
    lc_open : iterable
        Pseudo-observed light curves.

    Returns
    -------
    tuple
        (dt, t0) where dt is a list of lists (one per filter) and t0 contains
        the earliest detection time across filters.
    """
    # pylint: disable=too-many-locals,too-many-branches,too-many-statements,too-many-nested-blocks

    filter_color = ["b", "c", "g", "orange", "r", "m"]
    dt = [[] for _ in range(6)]
    t0 = []
    for lc in lc_open:
        if lc == 0:
            dt = [x.append(np.nan) or x for x in dt]
        else:
            mag_list = []
            time_list = []
            for f in filter_color:
                mag_list.append(
                    [
                        lc["mags"][j]
                        for j in range(len(lc["mags"]))
                        if lc["mags"][j] < lc["mags_lim"][j] and lc["filt"][j] == f
                    ]
                )
                time_list.append(
                    [
                        lc["time"][j]
                        for j in range(len(lc["mags"]))
                        if lc["mags"][j] < lc["mags_lim"][j] and lc["filt"][j] == f
                    ]
                )
            all_t0 = []
            for j in range(6):
                mag = mag_list[j]
                time = time_list[j]
                if mag:
                    all_t0.append(time[0])
                    dt[j].append(time[mag.index(min(mag))] - time[0])
                else:
                    dt[j].append(np.nan)
                    all_t0.append(np.nan)
            t0.append(np.nanmin(all_t0) if np.nanmin(all_t0) != np.nan else np.nan)
    return dt, t0


# pylint: disable=too-many-locals,too-many-branches,too-many-statements,too-many-nested-blocks
def rate(lc_open, data="simu"):
    """Compute decreasing and increasing rates for each pseudo-observed light curve.

    Returns three lists (per filter) containing the increasing rate, the rate in the
    first third of the decreasing part, and the rate in the final third.
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
            if data == "simu":
                filter_color = ["b", "c", "g", "orange", "r", "m"]
                for f in filter_color:
                    mag_list.append(
                        [
                            lc["mags"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"][j] and lc["filt"][j] == f
                        ]
                    )
                    time_list.append(
                        [
                            lc["time"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"][j] and lc["filt"][j] == f
                        ]
                    )
            elif data == "elasticc":
                filter_color = ["u", "g", "r", "i", "z", "Y"]
                for f in filter_color:
                    mag_list.append(
                        [
                            lc["mags"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"] and lc["filt"][j] == f
                        ]
                    )
                    time_list.append(
                        [
                            lc["time"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"] and lc["filt"][j] == f
                        ]
                    )
            for j in range(6):
                mag = mag_list[j]
                time = time_list[j]
                if len(mag) > 3 and mag.index(min(mag)) != len(mag) - 1:
                    mdec = mag[mag.index(min(mag)) : len(mag) - 1]
                    tdec = time[mag.index(min(mag)) : len(mag) - 1]
                    if mag.index(min(mag)) != 0:
                        minc = mag[0 : mag.index(min(mag))]
                        tinc = time[0 : mag.index(min(mag))]
                        dtinc = (
                            tinc[minc.index(min(minc))]
                            - tinc[minc.index(max(minc))]
                        )
                        rate_inc[j].append((min(minc) - max(minc)) / dtinc if dtinc > 1 else np.nan)
                    else:
                        rate_inc[j].append(np.nan)
                    if len(mdec) > 3:
                        m1 = mdec[0 : int(len(mdec) / 3)]
                        m3 = mdec[2 * int(len(mdec) / 3) : len(mdec)]
                        t1 = tdec[0 : int(len(tdec) / 3)]
                        t3 = tdec[2 * int(len(tdec) / 3) : len(tdec)]
                        dt1 = t1[m1.index(max(m1))] - t1[m1.index(min(m1))]
                        dt3 = t3[m3.index(max(m3))] - t3[m3.index(min(m3))]
                        if dt1 > 1 and dt3 > 1:
                            rate_dec_1[j].append((max(m1) - min(m1)) / dt1)
                            rate_dec_3[j].append((max(m3) - min(m3)) / dt3)
                        else:
                            rate_dec_1[j].append(np.nan)
                            rate_dec_3[j].append(np.nan)
                    else:
                        rate_dec_1[j].append(np.nan)
                        rate_dec_3[j].append(np.nan)
                elif len(mag) > 3 and mag.index(min(mag)) == len(mag) - 1:
                    minc = mag[0 : mag.index(min(mag))]
                    tinc = time[0 : mag.index(min(mag))]
                    dtinc = tinc[minc.index(min(minc))] - tinc[minc.index(max(minc))]
                    rate_inc[j].append((min(minc) - max(minc)) / dtinc if dtinc > 1 else np.nan)
                    rate_dec_1[j].append(np.nan)
                    rate_dec_3[j].append(np.nan)
                else:
                    rate_inc[j].append(np.nan)
                    rate_dec_1[j].append(np.nan)
                    rate_dec_3[j].append(np.nan)
    return rate_inc, rate_dec_1, rate_dec_3


def color(lc_open, data="simu"):
    """Calculate the mean g-r color for each pseudo-observed light curve.

    Returns a list of mean g-r values (``np.nan`` when not computable).
    """

    colors = []
    for lc in lc_open:
        if lc == 0:
            colors.append(np.nan)
        else:
            mag_list = []
            time_list = []
            if data == "simu":
                filter_color = ["c", "g"]
                for f in filter_color:
                    mag_list.append(
                        [
                            lc["mags"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"][j] and lc["filt"][j] == f
                        ]
                    )
                    time_list.append(
                        [
                            lc["time"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"][j] and lc["filt"][j] == f
                        ]
                    )
            elif data == "elasticc":
                filter_color = ["g", "r"]
                for f in filter_color:
                    mag_list.append(
                        [
                            lc["mags"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"] and lc["filt"][j] == f
                        ]
                    )
                    time_list.append(
                        [
                            lc["time"][j]
                            for j in range(len(lc["mags"]))
                            if lc["mags"][j] < lc["mags_lim"] and lc["filt"][j] == f
                        ]
                    )
            mag_r = mag_list[1]
            mag_g = mag_list[0]
            if (
                len(mag_r) != 0
                and len(mag_g) != 0
                and (mag_r.index(min(mag_r)) != len(mag_r) - 1)
                and (mag_g.index(min(mag_g)) != len(mag_g) - 1)
            ):
                time_r = time_list[1]
                time_g = time_list[0]
                mdec_r = mag_r[mag_r.index(min(mag_r)) : len(mag_r) - 1]
                tdec_r = time_r[mag_r.index(min(mag_r)) : len(mag_r) - 1]
                mdec_g = mag_g[mag_g.index(min(mag_g)) : len(mag_g) - 1]
                tdec_g = time_g[mag_g.index(min(mag_g)) : len(mag_g) - 1]
                if (
                    len(mdec_r) > 1
                    and (tdec_g[-1] > tdec_r[0])
                    and (tdec_r[-1] > tdec_g[0])
                    and (tdec_r[-1] - tdec_r[0] > tdec_g[-1] - tdec_g[0])
                ):
                    interpolation = scipy.interpolate.interp1d(tdec_r, mdec_r)
                    while tdec_g[0] < tdec_r[0]:
                        tdec_g.pop(0)
                        mdec_g.pop(0)
                    while tdec_g[-1] > tdec_r[-1]:
                        tdec_g.pop()
                        mdec_g.pop()
                    magr_int = interpolation(tdec_g)
                    colors.append(np.mean(mdec_g - magr_int))
                elif (
                    len(mdec_g) > 1
                    and (tdec_r[-1] > tdec_g[0])
                    and (tdec_g[-1] > tdec_r[0])
                    and (tdec_r[-1] - tdec_r[0] < tdec_g[-1] - tdec_g[0])
                ):
                    interpolation = scipy.interpolate.interp1d(tdec_g, mdec_g)
                    while tdec_r[0] < tdec_g[0]:
                        tdec_r.pop(0)
                        mdec_r.pop(0)
                    while tdec_r[-1] > tdec_g[-1]:
                        tdec_r.pop()
                        mdec_r.pop()
                    magg_int = interpolation(tdec_r)
                    colors.append(np.mean(magg_int - mdec_r))
                else:
                    colors.append(np.nan)
            else:
                colors.append(np.nan)
    return colors


# pylint: disable=too-many-arguments,too-many-positional-arguments

def heatmap(
    lc_open,
    mag_min,
    time_max,
    t0,
    rate_dec_1,
    rate_dec_3,
    rate_inc,
    color_vals,
    dt,
    data="simu",
    parameters="all",
    annot=True,
):
    """Calculate correlations and plot a heatmap between model parameters and features.

    Parameters
    ----------
    lc_open, mag_min, time_max, t0, rate_dec_1, rate_dec_3, rate_inc, color_vals, dt
        Pre-computed feature lists.
    data : str, optional
        ``"simu"`` or ``"elasticc"``.
    parameters : str, optional
        ``"model"``, ``"all"`` (default) or ``"features"``.
    annot : bool, optional
        Whether to annotate heatmap cells.
    """

    dt_mean = []
    rate_dec_1_mean = []
    rate_dec_3_mean = []
    rate_inc_mean = []
    all_results_not_zeros = None
    for i in range(len(lc_open)):
        if data == "simu":
            dt_one = [dt[j][i] for j in range(len(dt)) if not np.isnan(dt[j][i])]
            dt_mean.append(np.mean(dt_one))
            rate_dec_1_one = [
                rate_dec_1[j][i]
                for j in range(len(rate_dec_1))
                if not math.isnan(rate_dec_1[j][i])
            ]
            rate_dec_3_one = [
                rate_dec_3[j][i]
                for j in range(len(rate_dec_3))
                if not math.isnan(rate_dec_3[j][i])
            ]
            rate_inc_one = [
                rate_inc[j][i]
                for j in range(len(rate_inc))
                if not math.isnan(rate_inc[j][i])
            ]
            rate_dec_1_mean.append(np.mean(rate_dec_1_one))
            rate_dec_3_mean.append(np.mean(rate_dec_3_one))
            rate_inc_mean.append(np.mean(rate_inc_one))
        elif data == "elasticc":
            rate_dec_1_one = [
                rate_dec_1[j][i]
                for j in range(len(rate_dec_1))
                if not math.isnan(rate_dec_1[j][i])
            ]
            rate_dec_3_one = [
                rate_dec_3[j][i]
                for j in range(len(rate_dec_3))
                if not math.isnan(rate_dec_3[j][i])
            ]
            rate_inc_one = [
                rate_inc[j][i]
                for j in range(len(rate_inc))
                if not math.isnan(rate_inc[j][i])
            ]
            rate_dec_1_mean.append(np.mean(rate_dec_1_one))
            rate_dec_3_mean.append(np.mean(rate_dec_3_one))
            rate_inc_mean.append(np.mean(rate_inc_one))

    results = pd.DataFrame(
        list(
            zip(
                mag_min,
                time_max,
                t0,
                rate_dec_1_mean,
                rate_dec_3_mean,
                rate_inc_mean,
                color_vals,
                dt_mean,
            )
        ),
        columns=[
            "mag_peak",
            "time_peak",
            "first_detect",
            "rate_dec_1",
            "rate_dec_3",
            "rate_inc",
            "color",
            "dt",
        ],
    )

    results_not_zeros = results[(results["time_peak"] > 0) & (results["dt"] > 0)]
    correlations = None
    if data == "simu":
        lc_param = []
        dict_nan = {
            "jetType": -1,
            "specType": 0,
            "b": 4,
            "thetaObs": np.nan,
            "E0": np.nan,
            "thetaWing": np.nan,
            "thetaCore": np.nan,
            "n0": np.nan,
            "p": 2.2,
            "epsilon_e": 0.1,
            "epsilon_B": 0.01,
            "xi_N": 1.0,
            "d_L": 1.0e28,
            "z": np.nan,
        }
        for i in range(len(lc_open)):
            if lc_open[i] == 0:
                lc_param.append(dict_nan)
            else:
                lc_param.append(lc_open[i]["config"])
        df = pd.DataFrame(
            lc_param, columns=["E0", "n0", "z", "thetaObs", "thetaCore", "thetaWing"]
        )
        sns.set_theme(style="whitegrid")
        df["log10(E0)"] = np.log10(df["E0"])
        df["log10(n0)"] = np.log10(df["n0"])
        df["cos(thetaObs)"] = np.cos(df["thetaObs"])
        all_results = pd.concat([df, results], axis=1)
        display(all_results.dropna(axis=0))
        all_results_not_zeros = pd.concat([df, results_not_zeros], axis=1)
    elif data == "elasticc":
        all_results_not_zeros = results

    plt.figure(figsize=(15, 13))
    if parameters == "model":
        correlations = df.corr(numeric_only=True)
        sns.set(font_scale=1.7)
        sns.heatmap(
            df.corr(numeric_only=True), annot=annot, annot_kws={"size": 18}, center=0
        )
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        g = sns.PairGrid(df.corr(numeric_only=True), diag_sharey=False, corner=True)
        g.map_upper(sns.histplot)
        g.map_lower(sns.kdeplot, fill=True)
        g.map_diag(sns.kdeplot, color="midnightblue", lw=1.5)
        g.map_diag(sns.histplot, kde=True)
    elif parameters == "all":
        correlations = all_results_not_zeros.corr(numeric_only=True)
        sns.set(font_scale=1.3)
        sns.heatmap(
            all_results_not_zeros.corr(numeric_only=True),
            annot=annot,
            annot_kws={"size": 11},
            center=0,
        )
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        g = sns.PairGrid(
            all_results_not_zeros.corr(numeric_only=True),
            diag_sharey=False,
            corner=True,
        )
        g.map_lower(sns.kdeplot, fill=True)
        g.map_diag(sns.kdeplot, lw=1.5)
    elif parameters == "features":
        correlations = results_not_zeros.corr(numeric_only=True)
        sns.set(font_scale=1.7)
        sns.heatmap(
            results_not_zeros.corr(numeric_only=True),
            annot=annot,
            annot_kws={"size": 18},
            center=0,
        )
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        g = sns.PairGrid(
            results_not_zeros.corr(numeric_only=True), diag_sharey=False, corner=True
        )
        g.map_lower(sns.kdeplot, fill=True)
        g.map_diag(sns.kdeplot, lw=1.5)
    return correlations

# ---------------------------------------------------------------------------
# Deprecated CamelCase aliases – simple name bindings for backward compatibility.
# pylint: disable=invalid-name
# ---------------------------------------------------------------------------

MinimalMagnitude = minimal_magnitude
PeakTime = peak_time
DurationBetweenFirstAndPeak = duration_between_first_and_peak
Rate = rate
Color = color
Heatmap = heatmap

# End of file
