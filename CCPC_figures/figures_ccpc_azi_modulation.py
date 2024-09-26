import os
import glob
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
import seaborn as sns
import time
import logging
import python_scripts.find_mean_inc_angle
from python_scripts.find_mean_inc_angle import find_mean_incidence_angle
from python_scripts.wind_inversion_functions import cmod5n_forward

### To make maps with Cartopy using local data
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader

# TODO: introduce the future package name below
from plot_MACS_modulation_azimu_wrt_lambda_windpseed import (
    load_and_merge_dataframes,
    add_ancillary_variables,
    mean_curve_calcandplot,
    mean_curve_calc2,
    mean_curve_calc_180_180
)

# Define local data path
# local_shapefile_dir = '/home1/datahome/ljessel/Cartopy/'  # Remplacez par le chemin vers vos fichiers décompressés
local_shapefile_dir = "/home/datawork-cersat-public/cache/project/sarwave/tools/landmask/cartopy/shapefiles/natural_earth/physical"
shapefile_path = os.path.join(local_shapefile_dir, "ne_110m_coastline.shp")

# Read the local data by using Reader
coastline = Reader(shapefile_path).geometries()
coastline_feature = cfeature.ShapelyFeature(coastline, ccrs.PlateCarree())

sat_colors = {"S1A": "blue", "S1B": "red"}
mean_iangle = np.array(
    [
        31.64091048,
        32.83585645,
        34.03080242,
        35.15106426,
        36.79411497,
        37.83969269,
        38.88527042,
        39.93084814,
        40.97642586,
        42.17137183,
        43.14226543,
        44.03847491,
        44.93468438,
    ]
)  # hard coded, since Alexis Lucas did somthing over complicated to find the incidence angles
mean_iangle_iw1 = mean_iangle[:4]
mean_iangle_iw2 = mean_iangle[4:9]
mean_iangle_iw3 = mean_iangle[9:]


def ccpc_single_verification_figure_azi_modulation(df, satellite, burstkind):
    """
    sanity check figure

    :param df: pandas.DataFrame containing data from .csv files
    :param satellite: str S1A or S1B or S1A+B
    :param burstkind: str intraburst or interburst
    :return:
    """
    ### CCPC ###
    mask_verif = (
        (df["Wspeed"].values > 13)
        & (df["Wspeed"].values < 17)
        & (df["incidence"].values > 39)
        & (df["incidence"].values < 41)
    )  # selection of a range of 2 m/s around 15m/s wind speed and a range of 1° around 40° of incidence
    # print("mask", mask_verif.size, mask_verif.sum())
    NRCS_verif = df["sigma0_dB_filt"][mask_verif]
    if burstkind == "intraburst":
        ccpc_Re_verif = df["CCPC_filt_Re"][mask_verif]
    elif burstkind == "interburst":
        ccpc_Re_verif = df["CCPC_overlap_filt_Re"][mask_verif]
    else:
        raise Exception("burstkind : %s not handled", burstkind)
    az_wdir_verif = (
        df["wdir_az"].loc[ccpc_Re_verif.index].values
    )  # selection of the az_wdir corresponding to the range
    # print(az_wdir_verif, az_wdir_verif.shape)
    # print(ccpc_Re_verif.values, ccpc_Re_verif.values.shape)
    # plot figure
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        left=0.1,
        right=0.9,
        bottom=0.1,
        top=0.9,
        wspace=0.05,
        hspace=0.05,
    )

    ymax = 0.08
    # xmin = 0
    # xmax = 360
    xmin = -5
    xmax = 365

    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)
    if satellite == 'S1A+B':
        color = 'green'
    else:
        color = sat_colors[satellite]
    ## Density Plot
    sns.kdeplot(
        x=az_wdir_verif,
        y=ccpc_Re_verif.values,
        color=color,
        fill=True,
        ax=ax,
        label=satellite,
        bw_adjust = .5
    )  #

    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel("Azimuthal wind direction [°]")
    ax.set_ylabel("CCPC -Real part - vv - [$m^{-6}$]")

    ax.hlines(0, -70, 420, color="black", lw=1)
    ax.vlines(
        [90, 270],
        -ymax,
        ymax,
        color="teal",
        label="crosswind",
        linestyles="dashdot",
        alpha=0.6,
        lw=2.5,
    )
    ax.vlines(
        180,
        -ymax,
        ymax,
        color="black",
        label="downwind",
        linestyles="dashdot",
        alpha=0.6,
        lw=2.5,
    )
    ax.vlines(
        [1, 359.9],
        -ymax,
        ymax,
        color="maroon",
        label="upwind",
        linestyles="dashdot",
        alpha=0.6,
        lw=3,
    )
    ax.hlines(0, -70, 420, color="black", lw=1)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-ymax, ymax)
    ax.grid(ls="--")

    ## Histogramms
    binwidth = 1
    x_bins = np.arange(0, 360 + binwidth, binwidth * 5)
    ax_histx.hist(az_wdir_verif, bins=x_bins, color=color)
    ax_histx.grid(ls="--")
    y_bins = np.arange(-ymax, ymax, binwidth * ymax ** 2)
    ax_histy.hist(ccpc_Re_verif, bins=y_bins, color=color, orientation="horizontal")
    ax_histy.grid(ls="--")

    ## box title
    txt_str = (
        "Number of points: %d \nwind : 15 $\pm$ 2 m/s \nincidence : 40 $\pm$ 1 ° \n%s %s (processing b07)"
        % (len(ccpc_Re_verif), satellite, burstkind)
    )
    props = dict(boxstyle="square", facecolor="white")
    ax.text(
        0.03,
        0.97,
        txt_str,
        transform=ax.transAxes,
        fontsize=8.5,
        verticalalignment="top",
        bbox=props,
    )

    # mean curve calculation
    mean_curve_calcandplot(
        az_wdir=az_wdir_verif,
        variable_tested=ccpc_Re_verif,
        lw=2,
        label="mean curve",
        ax=ax,
    )

    ax.legend()
    plt.show()


# def mean_curve_calc(az_wdir, variabletested):
#     # bin setting
#     delta_azi = 5
#     # bin_edges = np.arange(0, 361, delta_azi)  # Bins every 5 degrees
#     # bin_centers = bin_edges + delta_azi / 2  # Centre of each bin
#     # bin_centers = bin_centers[
#     #     :-1
#     # ]  # Delete the last element to keep the same length as bin_edges
#     bin_centers = np.arange(0,360,delta_azi)
#     valmean = np.zeros_like(bin_centers)  # Initialise mean Imacs values
#
#     for i in range(len(bin_centers)):
#         if i == 0 or i==len(bin_centers)-1: # special case to make to 0° and 360° be the same
#             mask = (az_wdir<=delta_azi) | (az_wdir>=360-delta_azi)
#         else:
#             mask = (az_wdir >= bin_centers[i] - delta_azi) & (
#                 az_wdir < bin_centers[i] + delta_azi
#             )  # Mask to selects the points in the bin
#         valmean[i] = (
#             np.mean(variabletested[mask]) if np.sum(mask) > 0 else np.nan
#         )  # Calculation of the Imacs mean value in the bin
#
#     return bin_centers, valmean

def convert_angles(angles):
    # Use modulo to ensure angles are within [0°, 360°] if needed.
    angles = np.mod(angles, 360)
    # Convert angles to [-180°, 180°] range
    angles = np.where(angles > 180, angles - 360, angles)
    return angles

def longepe_azi_figure(df, burstkind="intraburst"):
    """

    :param df: S1A+S1B
    :param burstkind:
    :return:
    """
    start_time = time.time()

    xmin = -180
    xmax = 180
    nb_pts = {}
    windpeeds = np.arange(2, 14, 2)
    windspeed_colors = ['b','g','r','c','m','y','k']
    incidences = [34.5,38.5,42.5]
    variables = {
        'Amplitude': {'ymin':0,'ymax':0.16}
        ,'Re': {'ymin':-0.15,'ymax':0.15}
        ,'Im': {'ymin':-0.07,'ymax':0.07}
    }
    # chosen_mean_iangle = mean_iangle_subswaths[subswath]
    fig, ax = plt.subplots(len(variables), len(incidences), figsize=(3*8, 3*6))
    if burstkind == "intraburst":

        varname_re = "CCPC_filt_Re"
        varname_im = "CCPC_filt_Im"
    elif burstkind == "interburst":
        varname_re = "CCPC_overlap_filt_Re"
        varname_im = "CCPC_overlap_filt_Im"
    values_ccpc = {'Amplitude': abs(df[varname_re]+1j*df[varname_im]), # TO BE CHECKED wrt paper
                   'Re':df[varname_re],
                   'Im':df[varname_im],}
    delta_ws = 2 #m/s
    delta_inc = 1 # degree
    for i,var_x in enumerate(variables):

        for j in range(len(incidences)):  # "loop over the incidence angles
            for wwi,wscenter in enumerate(windpeeds):
                wnd_mask = (df["Wspeed"] > wscenter-delta_ws) & (df["Wspeed"] < wscenter+delta_ws)
                inc_mask = (df["incidence"] > incidences[j]-delta_inc) & (df["incidence"] < incidences[j]+delta_inc)
                descending = df["ground_heading"][abs(df["ground_heading"]) > 150]
                ascending = df["ground_heading"][abs(df["ground_heading"]) < 30]
                # Selection of data
                ccpc_selection = values_ccpc[var_x][wnd_mask & inc_mask].dropna()
                # Imacs_sel = Imacs_sel[Imacs_sel != 0]  # selection of non nul values
                ccpc_sel_asc = ccpc_selection.loc[
                    ccpc_selection.index.intersection(ascending.index)
                ]
                ccpc_sel_desc = ccpc_selection.loc[
                    ccpc_selection.index.intersection(descending.index)
                ]
                nb_pts["_ascen"] = len(ccpc_sel_asc)
                nb_pts["_descen"] = len(ccpc_sel_desc)
                logging.debug('nb points asc : %s',nb_pts["_ascen"])
                az_wdir_sel_asc = df["wdir_az"].loc[
                    ccpc_sel_asc.index]
                az_wdir_sel_desc = df["wdir_az"].loc[
                    ccpc_sel_desc.index]
                az_wdir_sel_asc = convert_angles(az_wdir_sel_asc)
                az_wdir_sel_desc = convert_angles(az_wdir_sel_desc)
                assert (az_wdir_sel_asc>180).sum()==0
                # ## Dataframe boundaries extension for ascending data
                # # Create a new dataframe to duplicate data from the left (300° --> 0°)
                # Imacs_before = Imacs_sel_asc[az_wdir_sel_asc >= 300]
                # az_wdir_before = az_wdir_sel_asc[az_wdir_sel_asc >= 300]
                # az_wdir_before -= 360  # modify the azimutal wind direction to get into the [-60 - 0] range
                #
                # # Create a new dataframe to duplicate data from the right (360 --> 60)
                # Imacs_after = Imacs_sel_asc[az_wdir_sel_asc <= 60]
                # az_wdir_after = az_wdir_sel_asc[az_wdir_sel_asc <= 60]
                # az_wdir_after += 360  # modify the azimutal wind direction to get into the [360 - 60] range
                #
                # # Add the two duplicated datas to the original dataframe
                # imacs_ext_s1a_asc = pd.concat(
                #     [Imacs_before, Imacs_sel_asc, Imacs_after]
                # )
                # az_wdir_ext_s1a_asc = pd.concat(
                #     [az_wdir_before, az_wdir_sel_asc, az_wdir_after]
                # )
                #
                # ## Dataframe boundaries extension for descending data
                # # Create a new dataframe to duplicate data from the left (300° --> 0°)
                # Imacs_before = Imacs_sel_desc[az_wdir_sel_desc >= 300]
                # az_wdir_before = az_wdir_sel_desc[az_wdir_sel_desc >= 300]
                # az_wdir_before -= 360  # modify the azimutal wind direction to get into the [-60 - 0] range
                #
                # # Create a new dataframe to duplicate data from the right (360 --> 60)
                # Imacs_after = Imacs_sel_desc[az_wdir_sel_desc <= 60]
                # az_wdir_after = az_wdir_sel_desc[az_wdir_sel_desc <= 60]
                # az_wdir_after += 360  # modify the azimutal wind direction to get into the [360 - 60] range
                #
                # # Add the two duplicated datas to the original dataframe
                # imacs_ext_s1a_desc = pd.concat(
                #     [Imacs_before, Imacs_sel_desc, Imacs_after]
                # )
                # az_wdir_ext_s1a_desc = pd.concat(
                #     [az_wdir_before, az_wdir_sel_desc, az_wdir_after]
                # )

                ### mean curve calculation
                # bin_centers, Imacs_mean = mean_curve_calc2(
                #     az_wdir_ext_s1a_asc, imacs_ext_s1a_asc
                # )
                bin_centers, ccpc_mean,_ = mean_curve_calc_180_180(
                    az_wdir_sel_asc, ccpc_sel_asc
                )
                ax[i][j].plot(
                    bin_centers,
                    ccpc_mean,
                    color=windspeed_colors[wwi],
                    label="%1.1f $\pm$%i m/s asc"%(wscenter,delta_ws),
                    linestyle="-",
                    lw=2,
                )  # plot the mean curve
                bin_centers, ccpc_mean,_ = mean_curve_calc_180_180(
                    az_wdir_sel_desc, ccpc_sel_desc
                )
                ax[i][j].plot(
                    bin_centers,
                    ccpc_mean,
                    color=windspeed_colors[wwi],
                    label="%1.1f $\pm$%i m/s desc"%(wscenter,delta_ws),
                    linestyle="--",
                    lw=2,
                )  # plot the mean curve

            # ax[i][j].hlines(0, -70, 420, color="black", lw=1)
            # ax[i][j].vlines(
            #     [90, 270],
            #     -ymax,
            #     ymax,
            #     color="teal",
            #     linestyles="dashdot",
            #     alpha=0.6,
            #     lw=2.5,
            # )  # label='crosswind'
            # ax[i][j].vlines(
            #     180, -ymax, ymax, color="black", linestyles="dashdot", alpha=0.6, lw=2.5
            # )  # label='downwind',
            # ax[i][j].vlines(
            #     [1, 359.9],
            #     -ymax,
            #     ymax,
            #     color="maroon",
            #     linestyles="dashdot",
            #     alpha=0.6,
            #     lw=3,
            # )  # label='upwind'
            ax[i][j].set_xticks([-180, -90, 0, 90, 180])
            ax[i][j].tick_params(axis="x", labelsize=15)
            ax[i][j].tick_params(axis="y", labelsize=15)
            ax[i][j].set_xlim(xmin, xmax)
            ax[i][j].set_ylim(variables[var_x]['ymin'], variables[var_x]['ymax'])
            if i == len(incidences) - 1:
                ax[i][j].set_xlabel("Azimuthal wind direction [°]", fontsize=15)
            else:
                ax[i][j].set_xlabel("")
            if j == 0:
                ax[i][j].set_ylabel("CCPC %s - vv - []"%var_x, fontsize=15)
            else:
                ax[i][j].set_ylabel("")
            ax[i][j].grid(linestyle="--", color="gray", alpha=0.9)
            ax[i][j].legend(fontsize=8, loc="upper right")
            txt_str = (
                    "Ascending pts num: %d \nDescending pts num: %d \nincidence : %.1f $\pm$ %i°"
                    % (
                        nb_pts["_ascen"],
                        nb_pts["_descen"],
                        incidences[j],delta_inc,
                    )
            )
            props = dict(boxstyle="square", facecolor="white")
            ax[i][j].text(
                0.03,
                0.97,
                txt_str,
                transform=ax[i][j].transAxes,
                fontsize=12,
                verticalalignment="top",
                bbox=props,
            )
    fig.suptitle(
        "S1A+B CCPC versus azimuth wind direction | Processing B07 | %s"
        % (burstkind),
        y=0.93,
        fontsize=25,
    )

    end_time = time.time()
    print("Ploting time :", end_time - start_time, "s")

    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/intra/comp_S1AB/IMACS_comp_s1ab_ascdesc_iw1_inter.png')
    fig.show()





def asc_desc_ccpc_azi_modulation(
    dfs, part="Re", burstkind="intraburst", subswath="iw1"
):

    start_time = time.time()

    # ymin=-0.02;ymax=0.02 # interburst
    ymin = -0.08
    ymax = 0.08  # intraburst
    xmin = 0
    xmax = 360
    nb_pts = {}

    if burstkind == "intraburst":
        varname = "CCPC_filt_" + part
    elif burstkind == "interburst":
        varname = "CCPC_overlap_filt_" + part
    else:
        raise Exception("burstkind : %s not handled", burstkind)

    mean_iangle_subswaths = {
        "iw1": mean_iangle_iw1,
        "iw2": mean_iangle_iw2,
        "iw3": mean_iangle_iw3,
    }
    chosen_mean_iangle = mean_iangle_subswaths[subswath]
    fig, ax = plt.subplots(4, len(chosen_mean_iangle), figsize=(40, 25))
    for i in range(4):  # hard coded... not a good idea here
        for j in range(len(chosen_mean_iangle)):  #"loop over the incidence angle of a given subswath
            for sat in ["S1A", "S1B"]:

                df = dfs[sat]
                wnd_spd_5ms = df["Wspeed"][(df["Wspeed"] > 3) & (df["Wspeed"] < 7)]
                wnd_spd_10ms = df["Wspeed"][(df["Wspeed"] > 8) & (df["Wspeed"] < 12)]
                wnd_spd_15ms = df["Wspeed"][(df["Wspeed"] > 13) & (df["Wspeed"] < 17)]
                wnd_spd_20ms = df["Wspeed"][(df["Wspeed"] > 18) & (df["Wspeed"] < 22)]

                inc_tot = [
                    df["incidence"].loc[wnd_spd_5ms.index],
                    df["incidence"].loc[wnd_spd_10ms.index],
                    df["incidence"].loc[wnd_spd_15ms.index],
                    df["incidence"].loc[wnd_spd_20ms.index],
                ]
                # macs_Im_tot_s1a = [macs_Im_s1a.loc[wnd_spd_5ms.index],macs_Im_s1a.loc[wnd_spd_10ms.index],
                #     macs_Im_s1a.loc[wnd_spd_15ms.index], macs_Im_s1a.loc[wnd_spd_20ms.index]]
                if part == "real":
                    ccpc_re_tot_s1a = [
                        df[varname].loc[wnd_spd_5ms.index],
                        df[varname].loc[wnd_spd_10ms.index],
                        df[varname].loc[wnd_spd_15ms.index],
                        df[varname].loc[wnd_spd_20ms.index],
                    ]
                else:
                    ccpc_re_tot_s1a = [
                        df[varname].loc[wnd_spd_5ms.index],
                        df[varname].loc[wnd_spd_10ms.index],
                        df[varname].loc[wnd_spd_15ms.index],
                        df[varname].loc[wnd_spd_20ms.index],
                    ]
                descending = df["ground_heading"][abs(df["ground_heading"]) > 150]
                ascending = df["ground_heading"][abs(df["ground_heading"]) < 30]
                # Selection of data
                Imacs_sel = ccpc_re_tot_s1a[i][
                    (inc_tot[i] > chosen_mean_iangle[j] - 1)
                    & (inc_tot[i] < chosen_mean_iangle[j] + 1)
                ].dropna()  # selection of a range of 1° around one mean incidence angle
                Imacs_sel = Imacs_sel[Imacs_sel != 0]  # selection of non nul values
                Imacs_sel_asc = Imacs_sel.loc[
                    Imacs_sel.index.intersection(ascending.index)
                ]
                Imacs_sel_desc = Imacs_sel.loc[
                    Imacs_sel.index.intersection(descending.index)
                ]
                nb_pts[sat + "_ascen"] = len(Imacs_sel_asc)
                nb_pts[sat + "_descen"] = len(Imacs_sel_desc)
                az_wdir_sel_asc = df["wdir_az"].loc[
                    Imacs_sel_asc.index
                ]  # selection of the az_wdir corresponding to the range
                az_wdir_sel_desc = df["wdir_az"].loc[Imacs_sel_desc.index]

                ## Dataframe boundaries extension for ascending data
                # Create a new dataframe to duplicate data from the left (300° --> 0°)
                Imacs_before = Imacs_sel_asc[az_wdir_sel_asc >= 300]
                az_wdir_before = az_wdir_sel_asc[az_wdir_sel_asc >= 300]
                az_wdir_before -= 360  # modify the azimutal wind direction to get into the [-60 - 0] range

                # Create a new dataframe to duplicate data from the right (360 --> 60)
                Imacs_after = Imacs_sel_asc[az_wdir_sel_asc <= 60]
                az_wdir_after = az_wdir_sel_asc[az_wdir_sel_asc <= 60]
                az_wdir_after += 360  # modify the azimutal wind direction to get into the [360 - 60] range

                # Add the two duplicated datas to the original dataframe
                imacs_ext_s1a_asc = pd.concat(
                    [Imacs_before, Imacs_sel_asc, Imacs_after]
                )
                az_wdir_ext_s1a_asc = pd.concat(
                    [az_wdir_before, az_wdir_sel_asc, az_wdir_after]
                )

                ## Dataframe boundaries extension for descending data
                # Create a new dataframe to duplicate data from the left (300° --> 0°)
                Imacs_before = Imacs_sel_desc[az_wdir_sel_desc >= 300]
                az_wdir_before = az_wdir_sel_desc[az_wdir_sel_desc >= 300]
                az_wdir_before -= 360  # modify the azimutal wind direction to get into the [-60 - 0] range

                # Create a new dataframe to duplicate data from the right (360 --> 60)
                Imacs_after = Imacs_sel_desc[az_wdir_sel_desc <= 60]
                az_wdir_after = az_wdir_sel_desc[az_wdir_sel_desc <= 60]
                az_wdir_after += 360  # modify the azimutal wind direction to get into the [360 - 60] range

                # Add the two duplicated datas to the original dataframe
                imacs_ext_s1a_desc = pd.concat(
                    [Imacs_before, Imacs_sel_desc, Imacs_after]
                )
                az_wdir_ext_s1a_desc = pd.concat(
                    [az_wdir_before, az_wdir_sel_desc, az_wdir_after]
                )

                ### mean curve calculation
                bin_centers, Imacs_mean = mean_curve_calc2(
                    az_wdir_ext_s1a_asc, imacs_ext_s1a_asc
                )
                ax[i][j].plot(
                    bin_centers,
                    Imacs_mean,
                    color=sat_colors[sat],
                    label=sat + "$_{asc}$",
                    linestyle="-",
                    lw=2,
                )  # plot the mean curve
                bin_centers, Imacs_mean = mean_curve_calc2(
                    az_wdir_ext_s1a_desc, imacs_ext_s1a_desc
                )
                ax[i][j].plot(
                    bin_centers,
                    Imacs_mean,
                    color=sat_colors[sat],
                    label=sat + "$_{desc}$",
                    linestyle="--",
                    lw=2,
                )  # plot the mean curve

            ax[i][j].hlines(0, -70, 420, color="black", lw=1)
            ax[i][j].vlines(
                [90, 270],
                -ymax,
                ymax,
                color="teal",
                linestyles="dashdot",
                alpha=0.6,
                lw=2.5,
            )  # label='crosswind'
            ax[i][j].vlines(
                180, -ymax, ymax, color="black", linestyles="dashdot", alpha=0.6, lw=2.5
            )  # label='downwind',
            ax[i][j].vlines(
                [1, 359.9],
                -ymax,
                ymax,
                color="maroon",
                linestyles="dashdot",
                alpha=0.6,
                lw=3,
            )  #  label='upwind'
            ax[i][j].set_xticks([0, 90, 180, 270, 360])
            ax[i][j].tick_params(axis="x", labelsize=15)
            ax[i][j].tick_params(axis="y", labelsize=15)
            ax[i][j].set_xlim(xmin, xmax)
            ax[i][j].set_ylim(ymin, ymax)
            if i == len(inc_tot) - 1:
                ax[i][j].set_xlabel("Azimuthal wind direction [°]", fontsize=15)
            else:
                ax[i][j].set_xlabel("")
            if j == 0:
                ax[i][j].set_ylabel("CCPC Re - vv - []", fontsize=15)
            else:
                ax[i][j].set_ylabel("")
            ax[i][j].grid(linestyle="--", color="gray", alpha=0.9)
            ax[i][j].legend(fontsize=12, loc="upper right")
            txt_str = (
                "S1A$_{asc}$ pts num: %d \nS1A$_{desc}$ pts num: %d \nS1B$_{asc}$ pts num: %d \nS1B$_{desc}$ pts num: %d \nwind : %d $\pm$ 2 m/s \nincidence : %.1f $\pm$ 1° (IW1)"
                % (
                    nb_pts["S1A_ascen"],
                    nb_pts["S1A_descen"],
                    nb_pts["S1B_ascen"],
                    nb_pts["S1B_descen"],
                    (i + 1) * 5,
                    chosen_mean_iangle[j],
                )
            )
            props = dict(boxstyle="square", facecolor="white")
            ax[i][j].text(
                0.03,
                0.97,
                txt_str,
                transform=ax[i][j].transAxes,
                fontsize=12,
                verticalalignment="top",
                bbox=props,
            )
    fig.suptitle(
        "S1A & S1B comparison of mean CCPC %s part versus azimuth wind direction | Processing B07 \n%s subswath | %s"
        % (part, subswath.upper(), burstkind),
        y=0.93,
        fontsize=25,
    )

    end_time = time.time()
    print("Ploting time :", end_time - start_time, "s")

    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/intra/comp_S1AB/IMACS_comp_s1ab_ascdesc_iw1_inter.png')
    fig.show()
