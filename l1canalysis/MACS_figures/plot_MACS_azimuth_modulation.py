import time
from matplotlib import pyplot as plt
import pandas as pd
from l1canalysis.utils import (mean_iangle_iw1,mean_iangle_iw2,mean_iangle_iw3,sat_colors,mean_curve_calc2,
                               mean_curve_calcandplot,mean_and_std_curve_calc_180_180,subswath_colors,mean_iangle)
import logging
import numpy as np
import seaborn as sns

def verification_figure(df,satellite,product_id,ymax=0.1):
    ### IMACS ###
    mask_verif = (df['Wspeed'].values > 13) & (df['Wspeed'].values < 17) & (df['incidence'].values > 39) & (
                df['incidence'].values < 41)  # selection of a range of 2 m/s around 15m/s wind speed and a range of 1° around 40° of incidence
    logging.info('Nb points in mask_verif : %i',mask_verif.sum())
    # az_wdir = df['wdir_az']
    az_wdir = df['wdir_az_scat']
    logging.info('scat azimuth')
    macs_Im_verif = df['macs_Im_lambda_max=50.0'][mask_verif]
    az_wdir_verif = az_wdir.loc[macs_Im_verif.index]
    # plot figure
    fig = plt.figure(figsize=(8, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=(4, 1), height_ratios=(1, 4),
                          left=0.1, right=0.9, bottom=0.1, top=0.9,
                          wspace=0.05, hspace=0.05)

    #ymax = 0.03 # raw MACS
    #ymax = 0.1 # normalized

    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    ## Density Plot
    sns.kdeplot(x=az_wdir_verif, y=macs_Im_verif, color=sat_colors[satellite], fill=True, bw_adjust=.5, ax=ax, label=satellite)
    #ax.set_ylim(-0.03, 0.03)
    #ax.set_ylim(-0.1, 0.1)
    ax.set_xlabel('Azimuthal wind direction [°]')
    ax.set_ylabel('IMACS - vv - 50m [$m^{-6}$]')

    binwidth = 1  # for histogram side figures
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
        # mean curve calculation
        bin_centers, Imacs_mean,_,_ = mean_curve_calc2(az_wdir=az_wdir_verif, variabletested=macs_Im_verif)
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
        # mean curve calculation
        bin_centers, Imacs_mean, _, nb_values = mean_and_std_curve_calc_180_180(az_wdir=az_wdir_verif,
                                                                               variabletested=macs_Im_verif)
        azi0 = np.arange(-180, 181)
        x_bins = np.arange(-180, 180 + binwidth, binwidth * 5)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-ymax, ymax)
    ax.grid(ls='--')

    ## Histogramms
    ax_histx.hist(az_wdir_verif, bins=x_bins, color=sat_colors[satellite])
    ax_histx.grid(ls='--')
    y_bins = np.arange(-ymax, ymax, binwidth * ymax ** 2)
    ax_histy.hist(macs_Im_verif, bins=y_bins, color=sat_colors[satellite], orientation='horizontal')
    ax_histy.grid(ls='--')

    ## box title
    part2 = r'incidence : 40 $\pm$ 1 ° '
    part1 = r'wind : 15 $ \pm $ 2 m/s '
    part0 = 'Number of points: %d'%len(macs_Im_verif)
    part3 = satellite+' (processing %s)'%product_id
    txt_str = part0+'\n'+part1+'\n'+part2+'\n'+part3
    props = dict(boxstyle='square', facecolor='white')
    ax.text(0.03, 0.97, txt_str, transform=ax.transAxes, fontsize=8.5, verticalalignment='top', bbox=props)


    label = 'mean curve'
    ax.plot(bin_centers, Imacs_mean, label=label, linestyle='--', lw=2,
            color='black' )  # plot the mean curve

    ax.legend()
    plt.show()

def asc_desc_imacs_azi_modulation(
    dfs, part="Re", burstkind="intraburst", subswath="iw1",lambda_val='50',product_id='B09'
,ymax=0.4):

    start_time = time.time()

    #if burstkind == 'interburst':
    #    ymin = -0.02
    #    ymax = 0.02  # interburst
    #else:
    #    ymin = -0.04
    #    ymax = 0.04  # intraburst
    ymin = -ymax
    xmin = 0
    xmax = 360
    nb_pts = {}
    varname = 'macs_%s_lambda_max=%s'%(part,float(lambda_val))

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

                macs_s1 = [
                    df[varname].loc[wnd_spd_5ms.index],
                    df[varname].loc[wnd_spd_10ms.index],
                    df[varname].loc[wnd_spd_15ms.index],
                    df[varname].loc[wnd_spd_20ms.index],
                ]

                descending = df["ground_heading"][abs(df["ground_heading"]) > 150]
                ascending = df["ground_heading"][abs(df["ground_heading"]) < 30]
                # Selection of data
                Imacs_sel = macs_s1[i][
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
                varname_azimuth_angle = 'wdir_az_scat'
                az_wdir_sel_asc = df[varname_azimuth_angle].loc[
                    Imacs_sel_asc.index
                ]  # selection of the az_wdir corresponding to the range
                az_wdir_sel_desc = df[varname_azimuth_angle].loc[Imacs_sel_desc.index]

                ## Dataframe boundaries extension for ascending data
                # Create a new dataframe to duplicate data from the left (300° --> 0°)
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
                bin_centers, Imacs_mean, Imacs_std, nb_values = mean_curve_calc2(
                    az_wdir_sel_asc, Imacs_sel_asc
                )
                ax[i][j].plot(
                    bin_centers,
                    Imacs_mean,
                    color=sat_colors[sat],
                    label=sat + "$_{asc}$",
                    linestyle="-",
                    lw=2,
                )  # plot the mean curve
                bin_centers, Imacs_mean, Imacs_std, nb_values = mean_curve_calc2(
                    az_wdir_sel_desc, Imacs_sel_desc
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
                if part == 'Re':
                    var_name_label = 'MACS'
                else:
                    var_name_label = 'IMACS'
                ax[i][j].set_ylabel('%s %s - vv - %s [$m^{-6}$]' % (var_name_label, part, lambda_val), fontsize=15)
            else:
                ax[i][j].set_ylabel("")
            ax[i][j].grid(linestyle="--", color="gray", alpha=0.9)
            ax[i][j].legend(fontsize=12, loc="upper right")
            part1 = "S1A$_{asc}$ pts num: %d \nS1A$_{desc}$ pts num: %d \nS1B$_{asc}$ pts num: %d \nS1B$_{desc}$ pts num: %d"%(  nb_pts["S1A_ascen"],
                    nb_pts["S1A_descen"],
                    nb_pts["S1B_ascen"],
                    nb_pts["S1B_descen"],)
            part2 = r"wind : %d $\pm$ 2 m/s "%((i + 1) * 5)
            part3 = r"incidence : %.1f $\pm$ 1° %s"%(chosen_mean_iangle[j],subswath)
            txt_str = part1+'\n'+part2+'\n'+part3
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
        "S1A & S1B comparison of mean MACS %s part versus azimuth wind direction | Processing %s \n%s subswath | %s"
        % (part, product_id, subswath.upper(), burstkind),
        y=0.93,
        fontsize=25,
    )

    end_time = time.time()
    print("Ploting time :", end_time - start_time, "s")

    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/intra/comp_S1AB/IMACS_comp_s1ab_ascdesc_iw1_inter.png')
    fig.show()

def macs_azimutal_modulation_with_std(dfs,subswath="iw1",lambda_val='50',part="Re", burstkind="intraburst",ymax=0.4,product_id='B09'):
    """

    :param dfs: dict with keys S1A and S1B
    :param subswath: str iw1 iw2 or iw3
    :param lambda_val: str 50 up to 275
    :param part: Re or Im
    :param burstkind:
    :return:
    """
    start_time = time.time()
    #if burstkind =='interburst':
    #    ymin=-0.02
    #    ymax=0.02 # interburst
    #else:
        #ymin = -0.04
        #ymax = 0.04  # intraburst
    ymin = -ymax
    mean_iangle_subswaths = {
        "iw1": mean_iangle_iw1,
        "iw2": mean_iangle_iw2,
        "iw3": mean_iangle_iw3,
    }
    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))
    chosen_mean_iangle = mean_iangle_subswaths[subswath]
    fig, ax = plt.subplots(4, len(chosen_mean_iangle), figsize=(40, 25))
    for i in range(4):  # loop over the wind speed ranges. hard coded... not a good idea here
        for j in range(len(chosen_mean_iangle)):  # "loop over the incidence angle of a given subswath
            nb_pts = {}
            for sarunit in ['S1A','S1B']:
                df = dfs[sarunit]
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
                macs_s1 = [
                    df[varname].loc[wnd_spd_5ms.index],
                    df[varname].loc[wnd_spd_10ms.index],
                    df[varname].loc[wnd_spd_15ms.index],
                    df[varname].loc[wnd_spd_20ms.index],
                ]
                ### Selection of S1A data ###

                Imacs_sel = macs_s1[i][
                    (inc_tot[i] > chosen_mean_iangle[j] - 1)
                    & (inc_tot[i] < chosen_mean_iangle[j] + 1)]

                Imacs_sel = Imacs_sel[Imacs_sel != 0]  # selection of non nul values
                # az_wdir_sel = df["wdir_az"].loc[
                #     Imacs_sel.index
                # ]  # selection of the az_wdir corresponding to the range
                az_wdir_sel = df["wdir_az_scat"].loc[
                    Imacs_sel.index
                ]  # selection of the az_wdir corresponding to the range

                ### mean curve calculation
                label = sarunit
                nb_pts[sarunit] = len(az_wdir_sel)
                if np.max(az_wdir_sel) > 180:
                    bin_centers, Imacs_mean, Imacs_std, nb_values = mean_curve_calc2(az_wdir=az_wdir_sel,
                                                                                     variabletested=Imacs_sel)
                else:
                    bin_centers, Imacs_mean, Imacs_std, nb_values = mean_and_std_curve_calc_180_180(az_wdir=az_wdir_sel,
                                                                                                    variabletested=Imacs_sel)
                ax[i][j].plot(bin_centers, Imacs_mean, label=label, linestyle='-', lw=2,color=sat_colors[sarunit],)  # plot the mean curve
                # bin_centers, Imacs_mean = mean_curve_calc(az_wdir_ext_s1a, imacs_ext_s1a)
                # ax[i][j].plot(bin_centers, Imacs_mean, color='blue', label='S1A', linestyle='-',
                              # lw=2)  # plot the mean curve
                # bin_centers, Imacs_std = std_curve_calc(az_wdir_ext_s1a, imacs_ext_s1a)
                ax[i][j].fill_between(bin_centers,
                                      Imacs_mean - Imacs_std,
                                      Imacs_mean + Imacs_std,
                                      color=sat_colors[sarunit], alpha=0.3)

                # bin_centers, Imacs_mean = mean_curve_calc(az_wdir_ext_s1b, imacs_ext_s1b)
                # ax[i][j].plot(bin_centers, Imacs_mean, color='red', label='S1B', linestyle='-', lw=2)  # plot the mean curve
                # bin_centers, Imacs_std = std_curve_calc(az_wdir_ext_s1b, imacs_ext_s1b)
                # ax[i][j].fill_between(bin_centers,
                #                       Imacs_mean - Imacs_std,
                #                       Imacs_mean + Imacs_std,
                #                       color='red', alpha=0.3)
            if np.max(az_wdir_sel)>180:
                xmin,xmax = (0,360)
                ax[i][j].axhline(0, color='black', lw=1)
                ax[i][j].vlines([90,270], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6, lw=2.5)
                ax[i][j].vlines(180, -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6, lw=2.5)

                ax[i][j].vlines([1,359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
                ax[i][j].set_xticks([0,90,180,270,360])
            else:
                xmin, xmax = (-180, 180)
                ax[i][j].axhline(0, color='black', lw=1)
                ax[i][j].vlines([90, -90], -ymax, ymax, color='teal', label='crosswind', linestyles="dashdot", alpha=0.6,
                          lw=2.5)
                ax[i][j].vlines([-180, 180], -ymax, ymax, color='black', label='downwind', linestyles="dashdot", alpha=0.6,
                          lw=2.5)
                ax[i][j].vlines([0], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
                ax[i][j].hlines(0, -70, 420, color='black', lw=1)
                ax[i][j].set_xticks([-180, -90, 0, 90, 180])
            ax[i][j].tick_params(axis='x', labelsize=15)
            ax[i][j].tick_params(axis='y', labelsize=15)
            ax[i][j].set_xlim(xmin, xmax)
            ax[i][j].set_ylim(ymin, ymax)
            if i == len(inc_tot) - 1:
                ax[i][j].set_xlabel('Azimuthal wind direction [°]', fontsize=15)
            else:
                ax[i][j].set_xlabel('')
            if j == 0:
                if part=='Re':
                    var_name_label = 'MACS'
                else:
                    var_name_label = 'IMACS'
                ax[i][j].set_ylabel('%s %s - vv - %s [$m^{-6}$]'%(var_name_label,part,lambda_val), fontsize=15)
            else:
                ax[i][j].set_ylabel('')
            ax[i][j].grid(linestyle='--', color='gray', alpha=0.9)
            ax[i][j].legend(fontsize=12, loc='upper right')
            part1 = 'S1A pts num: %d \nS1B pts num: %d '%(nb_pts['S1A'], nb_pts['S1B'])
            part2 = r'wind : %d $ \pm $ 2 m/s '%((i + 1) * 5)
            part3 = r'incidence : %.1f $\pm$ 1° (%s)'%(chosen_mean_iangle[j],subswath)
            txt_str = part1+'\n'+part2+'\n'+part3
            props = dict(boxstyle='square', facecolor='white')
            ax[i][j].text(0.03, 0.97, txt_str, transform=ax[i][j].transAxes, fontsize=12, verticalalignment='top',
                          bbox=props)

    fig.suptitle(
        'S1A & S1B comparison of mean IMACS versus azimuth wind direction | Processing %s \n%s subswath | %s'%(product_id,subswath,burstkind),
        y=0.93, fontsize=25)
    end_time = time.time()
    print('Ploting time :', end_time - start_time, 's')

    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/intra/comp_S1AB/IMACS_comp_s1ab_IW1_intra.png')
    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/inter/comp_S1AB/IMACS_comp_s1ab_IW1_inter.png')
    fig.show()


def macs_per_subswath_per_windspeed(df,satellite='S1A+B',part='Re',burstkind='intraburst',polarization='vv',lambda_val='50',ymax=0.4):
    """

    section 4.3 in notebook S1A_S1B_IMACS_mean_analysis.ipynb

    :param df:
    :param satellite:
    :param part:
    :param burstkind:
    :param polarization:
    :param subswath:
    :return:
    """
    mean_iangle_subswaths = {
        "iw1": mean_iangle_iw1,
        "iw2": mean_iangle_iw2,
        "iw3": mean_iangle_iw3,
    }
    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))

    import matplotlib.gridspec as gridspec

    #ymin = -0.02
    #ymax = 0.02
    ymin = -ymax
    xmin = 0
    xmax = 360
    ls = ['solid', 'dotted', 'dashed', 'dashdot']

    # Individual figures dimension
    fig_width_per_plot = 7.5
    fig_height_per_plot = 6

    # Total width and height of the figure
    total_width = fig_width_per_plot * 5  # 5 columns
    total_height = fig_height_per_plot * 3  # 3 lines
    fig = plt.figure(figsize=(total_width, total_height))

    # gridspec for the 3 lines with variables coclumns; use of 10 columns to alow a semi-column offset
    gs = gridspec.GridSpec(3, 10, figure=fig)

    # Ajouter les subplots pour chaque ligne
    axes = []
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

    macs_s1 = [
        df[varname].loc[wnd_spd_5ms.index],
        df[varname].loc[wnd_spd_10ms.index],
        df[varname].loc[wnd_spd_15ms.index],
        df[varname].loc[wnd_spd_20ms.index],
    ]
    count_figure = 0
    for ssi,subswath in enumerate(['iw1','iw2','iw3']):
        chosen_mean_iangle = mean_iangle_subswaths[subswath]
        for iangle in chosen_mean_iangle:
            # print('chosen_mean_iangle',chosen_mean_iangle)
            ind_iangle = list(chosen_mean_iangle).index(iangle)
            c = 0
            if ssi in [0,2]: # we want 4 figures for iw1 and iw3 but 5 for iw2
                ax = fig.add_subplot(gs[ssi, 2 * ind_iangle + 1:2 * ind_iangle + 3])  # 1/2 columns offset creation
            else: #case midle line
                #2 * ind_iangle:2 * ind_iangle + 2
                ax = fig.add_subplot(gs[ssi, 2 * ind_iangle:2 * ind_iangle + 2])  # 1/2 columns offset creation
            axes.append(ax)
            for i in range(len(inc_tot)): # loop over the wind speed ranges

                ### Selection of S1 data ###
                Imacs_sel = macs_s1[i][(inc_tot[i] > iangle - 1) & (inc_tot[i] < iangle + 1)].dropna()  # selection of a range of 1° around one mean incidence angle
                Imacs_sel = Imacs_sel[Imacs_sel != 0]  # selection of non nul values
                # Imacs_sel_s1 = Imacs_sel
                c += len(Imacs_sel)
                az_wdir_sel = df["wdir_az_scat"].loc[
                    Imacs_sel.index
                ]  # selection of the az_wdir corresponding to the range

                # ## Dataframe boundaries extension
                # # Create a new dataframe to duplicate data from the left (300° --> 0°)
                # Imacs_before = Imacs_sel[az_wdir_sel >= 300]
                # az_wdir_before = az_wdir_sel[az_wdir_sel >= 300]
                # az_wdir_before -= 360  # modify the azimutal wind direction to get into the [-60 - 0] range
                #
                # # Create a new dataframe to duplicate data from the right (360 --> 60)
                # Imacs_after = Imacs_sel[az_wdir_sel <= 60]
                # az_wdir_after = az_wdir_sel[az_wdir_sel <= 60]
                # az_wdir_after += 360  # modify the azimutal wind direction to get into the [360 - 60] range
                #
                # # Add the two duplicated datas to the original dataframe
                # imacs_ext_s1 = pd.concat([Imacs_before, Imacs_sel, Imacs_after])
                # az_wdir_ext_s1 = pd.concat([az_wdir_before, az_wdir_sel, az_wdir_after])

                ### mean curve calculation
                # bin_centers, Imacs_mean,Imacs_std,nb_values = mean_and_std_curve_calc_180_180(az_wdir_sel, variabletested=Imacs_sel)
                bin_centers, Imacs_mean, Imacs_std, nb_values = mean_curve_calc2(az_wdir_sel,
                                                                                                variabletested=Imacs_sel)
                # bin_centers, Imacs_mean = mean_curve_calc(az_wdir_ext_s1, imacs_ext_s1)
                axes[count_figure].plot(bin_centers, Imacs_mean, label=r'wind speed: %.1f $\pm$ 2 m/s' % ((i + 1) * 5),
                                      color=subswath_colors[subswath], linestyle=ls[i],
                                      lw=2)  # plot the mean curve # alpha=1 - (ind_iangle * 0.2)
            if np.max(az_wdir_sel)>180:
                axes[count_figure].hlines(0, -70, 420, color='black', lw=1)
                axes[count_figure].vlines([90, 270], -ymax, ymax, color='teal', label='crosswind', linestyles='dashdot',
                                        alpha=0.6, lw=2.5)
                axes[count_figure].vlines(180, -ymax, ymax, color='black', label='downwind', linestyles='dashdot', alpha=0.6,
                                        lw=2.5)
                axes[count_figure].vlines([1, 359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot",
                                        alpha=0.6, lw=3)
                axes[count_figure].set_xticks([0, 90, 180, 270, 360])
            else:
                axes[count_figure].axhline(0, color='black', lw=1)
                axes[count_figure].vlines([-90, 90], -ymax, ymax, color='teal', label='crosswind', linestyles='dashdot',
                                        alpha=0.6, lw=2.5)
                axes[count_figure].vlines([-180,180], -ymax, ymax, color='black', label='downwind', linestyles='dashdot', alpha=0.6,
                                        lw=2.5)
                axes[count_figure].vlines([0], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot",
                                        alpha=0.6, lw=3)
                axes[count_figure].set_xticks([-180, -90,0, 90, 180])
            axes[count_figure].tick_params(axis='x', labelsize=15)
            axes[count_figure].tick_params(axis='y', labelsize=15)
            axes[count_figure].set_xlim(xmin, xmax)
            axes[count_figure].set_ylim(ymin, ymax)
            axes[count_figure].set_xlabel('Azimuthal wind direction [°]', fontsize=15)
            if ind_iangle == 0:
                axes[count_figure].set_ylabel(r'$\Im$(MACS) - %s - %sm'%(polarization,lambda_val), fontsize=15)
            else:
                axes[count_figure].set_ylabel('')
            axes[count_figure].grid(linestyle='--', color='gray', alpha=0.9)
            axes[count_figure].legend(fontsize=12, loc='lower center', ncols=2)
            axes[count_figure].set_title(r'incidence angle: %.1f $\pm$ 1° | %d points' % (iangle, c), fontsize=15)
            count_figure += 1

    # Ajuster les marges entre les subplots
    plt.subplots_adjust(left=0.05, right=0.95, top=0.90, bottom=0.05, hspace=0.3, wspace=0.5)

    # Ajouter un sous-titre pour chaque ligne
    fig.text(0.5, 0.94, 'IW1', ha='center', fontsize=16, fontweight='bold')
    fig.text(0.5, 0.62, 'IW2', ha='center', fontsize=16, fontweight='bold')
    fig.text(0.5, 0.31, 'IW3', ha='center', fontsize=16, fontweight='bold')

    # Ajouter un titre global
    fig.suptitle(satellite+' IMACS versus azimuthal wind direction | wind speed dependency | %s' % burstkind,
                 fontsize=20, fontweight='bold')
    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/inter/S1A+B/IMACSvsAZ_inc_filt_s1ab_allIW_inter.png')
    # Afficher la figure
    plt.show()


def macs_az_windspeed_inc_recap(df,satellite='S1A+B',part='Re',burstkind='intraburst',polarization='vv',lambda_val='50',ymax=0.4):
    fig, ax = plt.subplots(1, 4, figsize=(60, 12))
    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))
    #ymin = -0.04
    #ymax = 0.04
    ymin = -ymax
    xmin = 0
    xmax = 360
    mean_iangle_s1_sel = mean_iangle[::2]

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

    macs_s1 = [
        df[varname].loc[wnd_spd_5ms.index],
        df[varname].loc[wnd_spd_10ms.index],
        df[varname].loc[wnd_spd_15ms.index],
        df[varname].loc[wnd_spd_20ms.index],
    ]


    for i in range(4): # loop over the wind speed ranges
        for iangle in mean_iangle_s1_sel:
            ### Selection of S1 data ###
            Imacs_sel = macs_s1[i][(inc_tot[i] > iangle - 1) & (inc_tot[
                                i] < iangle + 1)].dropna()  # selection of a range of 1° around one mean incidence angle
            Imacs_sel = Imacs_sel[Imacs_sel != 0]  # selection of non nul values
            Imacs_sel_s1 = Imacs_sel
            # az_wdir_sel = az_wdir_s1_intra.loc[Imacs_sel.index]  # selection of the az_wdir corresponding to the range
            varname_azimuth_angle = 'wdir_az_scat'
            az_wdir_sel = df[varname_azimuth_angle].loc[
                Imacs_sel_s1.index
            ]

            ## Dataframe boundaries extension
            # # Create a new dataframe to duplicate data from the left (300° --> 0°)
            # Imacs_before = Imacs_sel[az_wdir_sel >= 300]
            # az_wdir_before = az_wdir_sel[az_wdir_sel >= 300]
            # az_wdir_before -= 360  # modify the azimutal wind direction to get into the [-60 - 0] range
            #
            # # Create a new dataframe to duplicate data from the right (360 --> 60)
            # Imacs_after = Imacs_sel[az_wdir_sel <= 60]
            # az_wdir_after = az_wdir_sel[az_wdir_sel <= 60]
            # az_wdir_after += 360  # modify the azimutal wind direction to get into the [360 - 60] range
            #
            # # Add the two duplicated datas to the original dataframe
            # imacs_ext_s1 = pd.concat([Imacs_before, Imacs_sel, Imacs_after])
            # az_wdir_ext_s1 = pd.concat([az_wdir_before, az_wdir_sel, az_wdir_after])

            ### mean curve calculation
            # bin_centers, Imacs_mean = mean_curve_calc(az_wdir_ext_s1, imacs_ext_s1)
            bin_centers, Imacs_mean, Imacs_std, nb_values = mean_curve_calc2(az_wdir_sel,
                                                                             variabletested=Imacs_sel_s1)
            ax[i].plot(bin_centers, Imacs_mean, label=r'inc: %.1f $\pm$ 1°' % (iangle), linestyle='-',
                       lw=2)  # plot the mean curve

        ax[i].hlines(0, -70, 420, color='black', lw=2)
        ax[i].vlines([90, 270], -ymax, ymax, color='teal', label='crosswind', linestyles='dashdot', alpha=0.6, lw=2.5)
        ax[i].vlines(180, -ymax, ymax, color='black', label='downwind', linestyles='dashdot', alpha=0.6, lw=2.5)
        ax[i].vlines([1, 359.9], -ymax, ymax, color='maroon', label='upwind', linestyles="dashdot", alpha=0.6, lw=3)
        ax[i].set_xticks([0, 90, 180, 270, 360])
        ax[i].tick_params(axis='x', labelsize=25)
        ax[i].tick_params(axis='y', labelsize=25)
        ax[i].set_xlim(xmin, xmax)
        ax[i].set_ylim(ymin, ymax)
        ax[i].set_xlabel('Azimuthal wind direction [°]', fontsize=25)
        if i == 0:
            ax[i].set_ylabel('IMACS - %s - %sm [$m^{-6}$]'%(polarization,lambda_val), fontsize=25)
        else:
            ax[i].set_ylabel('')
        ax[i].grid(linestyle='--', color='gray', alpha=0.9)
        ax[i].legend(fontsize=22, loc='lower center', ncols=3)
        ax[i].set_title(r'wind speed = %.1f $\pm$ 2 m/s | %d points' % ((i + 1) * 5, len(macs_s1[i])),
                        fontsize=25)
    fig.suptitle(satellite+' IMACS versus azimuthal wind direction | wind speed dependency \n %s %s'%(burstkind,polarization), fontsize=35)
    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/intra/S1A&B//IMACSvsAZ_wspd_filt_s1ab_allIW_intra.png')
    fig.show()
