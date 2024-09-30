### NRCS ###

import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from l1canalysis.miscellaneous.wind_inversion_functions import cmod5n_forward
from l1canalysis.MACS_figures.plot_MACS_modulation_azimu_wrt_lambda_windpseed import (
    mean_curve_calc2,
    mean_curve_calcandplot,
    mean_curve_calc_180_180
)
def nrcs_verification_figure(df,satellite):
    mask_verif = (df['Wspeed'].values > 13) & (df['Wspeed'].values < 17) & (df['Wspeed'].values > 39) & (
            df[
                'incidence'].values < 41)  # selection of a range of 2 m/s around 15m/s wind speed and a range of 1° around 40° of incidence
    NRCS_verif = df['sigma0_dB_filt'][mask_verif]
    az_wdir = df['wdir_az']
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
    xmin=0;xmax=360

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    ## density plots
    sns.kdeplot(x=az_wdir_verif, y=NRCS_verif, color='blue', fill=True, bw_adjust=.5, ax=ax, label='Sentinel-1A')
    ax.set_xlim(-5,365)
    ax.set_xlabel('Azimuthal wind direction [°]')
    ax.set_ylabel('NRCS [dB]')

    ax.hlines(0,-70,420, color='black', lw=1)
    ax.vlines([90,270], ymin, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines(180, ymin, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)
    ax.vlines([1,359.9], ymin, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
    ax.hlines(0,-70,420, color='black', lw=1)
    ax.set_xticks([0,90,180,270,360])
    ax.set_xlim(xmin,xmax); ax.set_ylim(ymin, ymax)
    ax.grid(ls='--')

    ## Histogramms
    binwidth = 1
    x_bins = np.arange(0,360 + binwidth, binwidth*5)
    ax_histx.hist(az_wdir_verif,bins=x_bins, color='blue')
    ax_histx.grid(ls='--')
    y_bins = np.arange(min(NRCS_verif), max(NRCS_verif) , binwidth/3)
    ax_histy.hist(NRCS_verif, bins=y_bins, color='blue', orientation='horizontal')
    ax_histy.grid(ls='--')


    ## box title
    txt_str = 'Number of points: %d \nwind : 15 $\pm$ 2 m/s \nincidence : 40 $\pm$ 1 ° \nSentinel-1A (processing b07)' % (len(NRCS_verif))
    props = dict(boxstyle = 'square', facecolor = 'white')
    ax.text(0.03, 0.97, txt_str, transform = ax.transAxes, fontsize = 8.5, verticalalignment = 'top', bbox = props)

    # mean curve calculation
    mean_curve_calcandplot(az_wdir_verif, NRCS_verif)

    # CMOD5 curve
    azi = np.arange(361)
    sig = cmod5n_forward(azi*0+15, azi, azi*0+40, True)
    ax.plot(azi,10*np.log10(sig), label='CMOD5',lw=3, color='firebrick')

    ax.legend()
    plt.show()