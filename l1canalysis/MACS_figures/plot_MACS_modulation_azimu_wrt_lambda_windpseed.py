"""
A Grouazel
creation: sept 2024
set of methods to analyze MACS as function of lambda max
"""

import numpy as np
import logging
from l1canalysis.utils import mean_and_std_curve_calc_180_180
from matplotlib import pyplot as plt






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
    # xmin=0
    # xmax=360
    xmin=-180
    xmax = 180
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
        bin_centers, valmean,valstd,nb_values = mean_and_std_curve_calc_180_180(az_wdir=az_wdir_verif,variabletested=macs_Im_verif)
        # mean_curve_calcandplot(ax=ax,az_wdir=az_wdir_verif, variable_tested=macs_Im_verif,label=label)
        ax.plot(bin_centers, valmean, label=label, linestyle='-', lw=2)  # plot the mean curve
    ax.set_xlim(-5,365)
    # ax.set_ylim(-0.05, 0.05)
    ax.set_xlabel('Azimuthal wind direction [°]')
    ax.set_ylabel('IMACS [$m^{-6}$]' )

    # ax.hlines(0,-70,420, color='black', lw=1)
    ax.axhline(0, color='black', lw=1)
    # ax.vlines([90,270], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([90, -90], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    # ax.vlines(180, -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([-180,180], -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    # ax.vlines([1,359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.vlines([0], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.hlines(0,-70,420, color='black', lw=1)
    # ax.set_xticks([0,90,180,270,360])
    ax.set_xticks([-180, -90, 0,90,180])
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
    # xmin=0
    # xmax=360
    xmin=-180
    xmax=180
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
        bin_centers, valmean, valstd, nb_values = mean_and_std_curve_calc_180_180(az_wdir=az_wdir_verif, variabletested=macs_Im_verif)
        ax.plot(bin_centers, valmean, label=label, linestyle='-', lw=2)  # plot the mean curve
        # mean_curve_calcandplot(ax=ax,az_wdir=az_wdir_verif, imacs=macs_Im_verif,label=label)
    # ax.set_xlim(-5,365)
    # ax.set_ylim(-0.05, 0.05)
    ax.set_xlabel('Azimuthal wind direction [°]')
    ax.set_ylabel('IMACS [$m^{-6}$]' )

    # ax.hlines(0,-70,420, color='black', lw=1)
    ax.axhline(0, color='black', lw=1)
    # ax.vlines([90,270], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([90, -90], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    # ax.vlines(180, -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([-180, 180], -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    # ax.vlines([1,359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.vlines([0], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.hlines(0, -70, 420, color='black', lw=1)
    # ax.set_xticks([0,90,180,270,360])
    ax.set_xticks([-180, -90, 0, 90, 180])
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








