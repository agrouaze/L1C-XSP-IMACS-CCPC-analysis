"""
A Grouazel
creation: sept 2024
set of methods to analyze MACS as function of lambda max
"""
import os
import glob
import pandas as pd
import time
import numpy as np
import logging
from aar.compute_angles import azimuth_direction,wind_direction_from_U_and_V
from matplotlib import pyplot as plt

def mean_curve_calcandplot(ax,az_wdir, variable_tested, label,lw=2):
    """

    :param ax: pyplot.axis
    :param az_wdir:np.ndarray
    :param variable_tested: np.ndarray imacs or ccpc
    :param label: str
    :param lw: int
    :return:
    """
    # new bin setting
    # delta_azi = 5
    # # bin_edges = np.arange(0, 361, delta_azi)  # Bins every 5 degrees
    # # bin_centers = bin_edges + delta_azi / 2   # Centre of each bin
    # # bin_centers = bin_centers[:-1]    # Del the last element to keep the same length as bin_edges
    # bin_centers = np.arange(0, 360, delta_azi)
    # mean_val = np.zeros_like(bin_centers)    # Initialise mean Imacs values
    #
    # for i in range(len(bin_centers)):
    #     # mask = (az_wdir >= bin_centers[i] - 5) & (az_wdir < bin_centers[i] + 5)     # Mask to selects the points in the bin
    #     if i == 0 or i==len(bin_centers)-1: # special case to make to 0° and 360° be the same
    #         mask = (az_wdir<=delta_azi) | (az_wdir>=360-delta_azi)
    #         logging.info('special bin #:%i bound az : %s %s',i,delta_azi,360-delta_azi)
    #     else:
    #         mask = (az_wdir >= bin_centers[i] - delta_azi) & (
    #             az_wdir < bin_centers[i] + delta_azi
    #         )  # Mask to selects the points in the bin
    #         logging.info('bound az : %s %s', bin_centers[i] - delta_azi, bin_centers[i] + delta_azi)
    #     mean_val[i] = np.mean(variable_tested[mask]) if np.sum(mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin
    bin_centers,mean_val = mean_curve_calc2(az_wdir, variabletested=variable_tested)
    ax.plot(bin_centers, mean_val, label=label, linestyle='-', lw = lw)     # plot the mean curve

def mean_curve_calc2(az_wdir,variabletested)->(np.ndarray,np.ndarray):
    """

    :param az_wdir: np.ndarray
    :param variabletested: np.ndarray
    :return:
    """
    delta_azi = 5
    bin_centers = np.arange(delta_azi / 2, 360 + delta_azi / 2, delta_azi)
    valmean = np.zeros_like(bin_centers)  # Initialise mean Imacs values
    for i in range(len(bin_centers)):
        first_bound = (bin_centers[i] - delta_azi) % 360
        last_bound = (bin_centers[i] + delta_azi) % 360
        # logging.info('first_bound %s  last_bound %s', first_bound,last_bound)
        if first_bound > last_bound:  # eg first=355 and last=5
            mask = (az_wdir >= first_bound) | (az_wdir < last_bound)
        else:  # eg first=25 and last=35
            mask = (az_wdir >= first_bound) & (az_wdir < last_bound)
        valmean[i] = (
            np.mean(variabletested[mask]) if np.sum(mask) > 0 else np.nan
        )  # Calculation of the Imacs mean value in the bin
    return bin_centers, valmean

def mean_curve_calc_180_180(az_wdir,variabletested)->(np.ndarray,np.ndarray,np.ndarray):
    """

    :param az_wdir: np.ndarray
    :param variabletested: np.ndarray
    :return:
    """
    delta_azi = 5
    bin_centers = np.arange(-180, 180, delta_azi)
    valmean = np.zeros_like(bin_centers).astype(np.float64)  # Initialise mean Imacs values
    nb_values = np.zeros_like(bin_centers).astype(np.float64)
    for i in range(len(bin_centers)):
        first_bound = (bin_centers[i] - delta_azi)
        if first_bound<-180:
            first_bound = first_bound%180
        last_bound = (bin_centers[i] + delta_azi)
        if last_bound>180:
            last_bound = -180+last_bound%180
        logging.debug('first_bound %s  last_bound %s', first_bound,last_bound)
        if first_bound > last_bound:  # eg first=170 and last=-165 or first=180 and last=-162
            mask = (az_wdir >= first_bound) | (az_wdir < last_bound)
        else:  # eg first=25 and last=35
            mask = (az_wdir >= first_bound) & (az_wdir < last_bound)
        if np.sum(mask) > 0:
            valmean[i] = np.mean(variabletested[mask])
        else:
            valmean[i] = np.nan
        nb_values[i] = mask.sum()
    return bin_centers, valmean,nb_values


def macs_az_modulation_lambda(df,satellite,burstkind,polarization,wspeedmax=17,wspeedmin=13,incmin=39,incmax=41,product_id='B07',subplot=None,fig=None,disp_legend=False):
    """

    :param df:
    :param satellite:
    :param burstkind:
    :param polarization:
    :param wspeedmax:
    :param wspeedmin:
    :param incmin:
    :param incmax:
    :param product_id:
    :return:
    """


    mask_verif = (df['Wspeed'].values > wspeedmin) & (df['Wspeed'].values < wspeedmax) & (df['incidence'].values > incmin) & (df['incidence'].values < incmax)  # selection of a range of 2 m/s around 15m/s wind speed and a range of 1° around 40° of incidence
    NRCS_verif = df['sigma0_dB_filt'][mask_verif]
    # macs_Im_verif = macs_Im_s1b[mask_verif]
    az_wdir_verif = df['wdir_az'].loc[NRCS_verif.index]       # selection of the az_wdir corresponding to the range
    ymax= 0.065
    xmin=0
    xmax=360
    if subplot is None:
        fig = plt.figure(figsize=(8, 6))
        # gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
        #                       left=0.1, right=0.9, bottom=0.1, top=0.9,
        #                       wspace=0.05, hspace=0.05)
        gs = fig.add_gridspec(1, 1)


        # Create the Axes.
        ax = fig.add_subplot(gs[0, 0])
    else:
        ax = fig.add_subplot(subplot)
    # ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    # ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)

    # no labels
    # ax_histx.tick_params(axis="x", labelbottom=False)
    # ax_histy.tick_params(axis="y", labelleft=False)

    ## Density Plot
    for lambda_val in np.arange(50,300,25):
        macs_Im_verif = df['macs_Im_lambda_max=%s'%float(lambda_val)][mask_verif]
        label = r'$\lambda_{max}$ = %s m'%(lambda_val)
        # sns.kdeplot(x=az_wdir_verif, y=macs_Im_verif, color='green', fill=True, bw_adjust=.5, ax=ax, label=label)

        # mean curve calculation
        mean_curve_calcandplot(ax=ax,az_wdir=az_wdir_verif, imacs=macs_Im_verif,label=label)
    ax.set_xlim(-5,365)
    # ax.set_ylim(-0.05, 0.05)
    ax.set_xlabel('Azimuthal wind direction [°]')
    ax.set_ylabel('IMACS [$m^{-6}$]' )

    ax.hlines(0,-70,420, color='black', lw=1)
    ax.vlines([90,270], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines(180, -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([1,359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.hlines(0,-70,420, color='black', lw=1)
    ax.set_xticks([0,90,180,270,360])
    ax.set_xlim(xmin,xmax); ax.set_ylim(-ymax, ymax)
    ax.grid(ls='--')

    ## Histogramms
    # binwidth = 1
    # x_bins = np.arange(0,360 + binwidth, binwidth*5)
    # ax_histx.hist(az_wdir_verif,bins=x_bins, color='green')
    # ax_histx.grid(ls='--')
    # y_bins = np.arange(-ymax, ymax , binwidth*ymax**2)
    # ax_histy.hist(macs_Im_verif,bins=y_bins, color='green', orientation='horizontal')
    # ax_histy.grid(ls='--')

    ## box title
    # txt_str = 'Number of points: %d \nwind : %s m/s to %s m/s \nincidence : 40 $\pm$ 1 ° \n%s (processing b07)' % (,satellite)
    if subplot is None:
        tit = '%s IW %s polarization: %s \n wind : %s m/s to %s m/s \nincidence : %s° to %s° \n product ID: %s nb points: %i'%(satellite,burstkind,polarization,wspeedmin,wspeedmax,incmin,incmax,product_id,mask_verif.sum())
        plt.title(tit)
    else:
        txt_str = 'wind : %s m/s to %s m/s \nincidence : %s° to %s° \n product ID: %s nb points: %i'%(wspeedmin,wspeedmax,incmin,incmax,product_id,mask_verif.sum())
        props = dict(boxstyle = 'square', facecolor = 'white')
        ax.text(0.03, 0.97, txt_str, transform = ax.transAxes, fontsize = 10, verticalalignment = 'bottom',horizontalalignment='center', bbox = props)


    if disp_legend is True:
        ax.legend(bbox_to_anchor=(1,1.05))
    # plt.show()
    return fig



def macs_az_modulation_lambda_sea_state(df,satellite,burstkind,polarization,wspeedmax=17,wspeedmin=13,hsmin=0,hsmax=4,
t0m1min=2,t0m1max=13,subplot=None,fig=None,disp_legend=False,axis_color=None):
    """

    :param df:
    :param satellite:
    :param burstkind:
    :param polarization:
    :param wspeedmax:
    :param wspeedmin:
    :param hsmin:
    :param hsmax:
    :return:
    """


    mask_verif = (df['Wspeed'].values > wspeedmin) & (df['Wspeed'].values < wspeedmax) & (df['hs'].values > hsmin) & (df['hs'].values < hsmax) &\
      (df['t0m1'].values > t0m1min) & (df['t0m1'].values < t0m1max)# selection of a range of 2 m/s around 15m/s wind speed and a range of 1° around 40° of incidence
    NRCS_verif = df['sigma0_dB_filt'][mask_verif]
    # macs_Im_verif = macs_Im_s1b[mask_verif]
    az_wdir_verif = df['wdir_az'].loc[NRCS_verif.index]       # selection of the az_wdir corresponding to the range
    ymax= 0.065
    xmin=0
    xmax=360
    if subplot is None:
        fig = plt.figure(figsize=(8, 6))
        # gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
        #                       left=0.1, right=0.9, bottom=0.1, top=0.9,
        #                       wspace=0.05, hspace=0.05)
        gs = fig.add_gridspec(1, 1)


        # Create the Axes.
        ax = fig.add_subplot(gs[0, 0])
    else:
        ax = fig.add_subplot(subplot)
    # ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    # ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)

    # no labels
    # ax_histx.tick_params(axis="x", labelbottom=False)
    # ax_histy.tick_params(axis="y", labelleft=False)

    ## Density Plot
    for lambda_val in np.arange(50,300,25):
        macs_Im_verif = df['macs_Im_lambda_max=%s'%float(lambda_val)][mask_verif]
        label = r'$\lambda_{max}$ = %s m'%(lambda_val)
        # sns.kdeplot(x=az_wdir_verif, y=macs_Im_verif, color='green', fill=True, bw_adjust=.5, ax=ax, label=label)

        # mean curve calculation
        mean_curve_calcandplot(ax=ax,az_wdir=az_wdir_verif, imacs=macs_Im_verif,label=label)
    ax.set_xlim(-5,365)
    # ax.set_ylim(-0.05, 0.05)
    ax.set_xlabel('Azimuthal wind direction [°]')
    ax.set_ylabel('IMACS [$m^{-6}$]' )

    ax.hlines(0,-70,420, color='black', lw=1)
    ax.vlines([90,270], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines(180, -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([1,359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.hlines(0,-70,420, color='black', lw=1)
    ax.set_xticks([0,90,180,270,360])
    ax.set_xlim(xmin,xmax); ax.set_ylim(-ymax, ymax)
    ax.grid(ls='--')
    if axis_color is not None:
        ax.tick_params(axis='x', colors=axis_color)
        ax.tick_params(axis='y', colors=axis_color)
        # ax.yaxis.label.set_color(axis_color)
        # ax.xaxis.label.set_color(axis_color)

    ## Histogramms
    # binwidth = 1
    # x_bins = np.arange(0,360 + binwidth, binwidth*5)
    # ax_histx.hist(az_wdir_verif,bins=x_bins, color='green')
    # ax_histx.grid(ls='--')
    # y_bins = np.arange(-ymax, ymax , binwidth*ymax**2)
    # ax_histy.hist(macs_Im_verif,bins=y_bins, color='green', orientation='horizontal')
    # ax_histy.grid(ls='--')

    ## box title
    # txt_str = 'Number of points: %d \nwind : %s m/s to %s m/s \nincidence : 40 $\pm$ 1 ° \n%s (processing b07)' % (,satellite)
    if subplot is None:
        tit = '%s IW %s polarization: %s \n wind : %s m/s to %s m/s \nHs : %sm to %sm Period : %s s to %s s \n nb points: %i'%(satellite,
        burstkind,polarization,wspeedmin,wspeedmax,hsmin,hsmax,t0m1min,t0m1max,mask_verif.sum())
        plt.title(tit)
    else:
        txt_str = 'wind : %s m/s to %s m/s \nHs : %sm to %sm Period : %s s to %s s\nnb points: %i'%(wspeedmin,
        wspeedmax,hsmin,hsmax,t0m1min,t0m1max,mask_verif.sum())
        props = dict(boxstyle = 'square', facecolor = 'white')
        ax.text(0.03, 0.97, txt_str, transform = ax.transAxes, fontsize = 10, verticalalignment = 'bottom',horizontalalignment='center', bbox = props)


    if disp_legend is True:
        ax.legend(bbox_to_anchor=(1,1.05))
    # plt.show()
    return fig








