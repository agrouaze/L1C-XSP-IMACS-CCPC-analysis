### NRCS ###

import numpy as np
import logging
from matplotlib import pyplot as plt
import seaborn as sns
from l1canalysis.miscellaneous.wind_inversion_functions import cmod5n_forward
from l1canalysis.utils import sat_colors,mean_and_std_curve_calc_180_180,mean_curve_calc2,mean_curve_calcandplot

def nrcs_verification_figure(df,satellite,azimuth_varname='wdir_az',density=True):
    """

    :param df: pd.DataFrame
    :param satellite: str S1A or S1B or S1A+B
    :return:
    """
    mask_verif = (df['Wspeed'].values > 13) & (df['Wspeed'].values < 17) & (df['incidence'].values > 39) & (
            df['incidence'].values < 41)  # selection of a range of 2 m/s around 15m/s wind speed and a range of 1° around 40° of incidence
    logging.info('mask_verif : %s',mask_verif.sum())
    NRCS_verif = df['sigma0_dB_filt'][mask_verif]
    # az_wdir = df['wdir_az']
    az_wdir = df[azimuth_varname]
    logging.info('azimuth_varname %s',azimuth_varname)
    # macs_Im_verif = df['macs_Im_lambda_max=50.0'][mask_verif]
    az_wdir_verif = az_wdir.loc[NRCS_verif.index]

    # plot figure
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(2, 2,  width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    # Create the Axes
    ax = fig.add_subplot(gs[1, 0])
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
    ymax=-7; ymin=-19


    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    if density is True:
        ## density plots
        sns.kdeplot(x=az_wdir_verif, y=NRCS_verif, color=sat_colors[satellite], fill=True, bw_adjust=.5, ax=ax, label=satellite)
    ax.set_xlabel('Azimuthal wind direction [°] \n %s'%azimuth_varname)
    ax.set_ylabel('NRCS [dB]')

    ax.axhline(0, color='black', lw=1)
    binwidth = 1 # for histogram side figures
    if np.max(az_wdir_verif) > 180:
        xmin = 0
        xmax = 360
        ax.hlines(0, -70, 420, color='black', lw=1)
        ax.vlines([90, 270], -ymax, ymax, color='teal', label='crosswind', linestyles='dashdot',
                                alpha=0.6, lw=2.5)
        ax.vlines(180, -ymax, ymax, color='black', label='downwind', linestyles='dashdot', alpha=0.6,
                                lw=2.5)
        ax.vlines([1, 359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot",
                                alpha=0.6, lw=3)
        ax.set_xticks([0, 90, 180, 270, 360])
        bin_centers, nrcs_mean,_,_ = mean_curve_calc2(az_wdir=az_wdir_verif, variabletested=NRCS_verif)
        azi0 = np.arange(361)
        x_bins = np.arange(0, 360 + binwidth, binwidth * 5)
    else:
        xmin = -180
        xmax = 180
        ax.axhline(0, color='black', lw=1)
        ax.vlines([-90, 90], -ymax, ymax, color='teal', label='crosswind', linestyles='dashdot',
                                alpha=0.6, lw=2.5)
        ax.vlines([-180, 180], -ymax, ymax, color='black', label='downwind', linestyles='dashdot',
                                alpha=0.6,
                                lw=2.5)
        ax.vlines([0], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot",
                                alpha=0.6, lw=3)
        ax.set_xticks([-180, -90, 0, 90, 180])
        bin_centers, nrcs_mean, _, nb_values = mean_and_std_curve_calc_180_180(az_wdir=az_wdir_verif,
                                                                               variabletested=NRCS_verif)
        azi0 = np.arange(-180, 181)
        x_bins = np.arange(-180, 180 + binwidth, binwidth * 5)
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin, ymax)
    ax.grid(ls='--')

    ## Histogramms

    ax_histx.hist(az_wdir_verif,bins=x_bins, color=sat_colors[satellite])
    ax_histx.grid(ls='--')
    y_bins = np.arange(min(NRCS_verif), max(NRCS_verif) , binwidth/3)
    ax_histy.hist(NRCS_verif, bins=y_bins, color=sat_colors[satellite], orientation='horizontal')
    ax_histy.grid(ls='--')


    ## box title
    txt_str = 'Number of points: %d \nwind : 15 $\pm$ 2 m/s \nincidence : 40 $\pm$ 1 ° \n%s (processing b07)' % (len(NRCS_verif),satellite)
    props = dict(boxstyle = 'square', facecolor = 'white')
    ax.text(0.03, 0.97, txt_str, transform = ax.transAxes, fontsize = 8.5, verticalalignment = 'top', bbox = props)

    # mean curve calculation
    # mean_curve_calcandplot(az_wdir_verif, NRCS_verif)

    label = 'mean'
    ax.plot(bin_centers, nrcs_mean, label=label, linestyle='-', lw=3,
            color='black', )  # plot the mean curve

    # CMOD5 curve
    # azi = np.arange(361)

    # azi = (azi0+360.)%360
    # import copy
    # azi = copy.copy(azi0)
    # azi[azi<0] = 360. - azi0[azi0<0]

    sig = cmod5n_forward(azi0*0+15, azi0, azi0*0+40, True)
    ax.plot(azi0,10*np.log10(sig), label='CMOD5',lw=3, color='firebrick')

    ax.legend()
    plt.show()