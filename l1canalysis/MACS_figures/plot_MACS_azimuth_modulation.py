import time
from matplotlib import pyplot as plt
import pandas as pd
from l1canalysis.utils import mean_iangle_iw1,mean_iangle_iw2,mean_iangle_iw3,sat_colors
from l1canalysis.MACS_figures.plot_MACS_modulation_azimu_wrt_lambda_windpseed import (
    mean_curve_calc2,
    mean_curve_calc_180_180
)
def asc_desc_imacs_azi_modulation(
    dfs, part="Re", burstkind="intraburst", subswath="iw1",lambda_val='50'
):

    start_time = time.time()

    # ymin=-0.02;ymax=0.02 # interburst
    ymin = -0.08
    ymax = 0.08  # intraburst
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
                ax[i][j].set_ylabel("MACS %s - vv - []"%part, fontsize=15)
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
        "S1A & S1B comparison of mean MACS %s part versus azimuth wind direction | Processing B07 \n%s subswath | %s"
        % (part, subswath.upper(), burstkind),
        y=0.93,
        fontsize=25,
    )

    end_time = time.time()
    print("Ploting time :", end_time - start_time, "s")

    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/IW_SLC_L1C_B07/intra/comp_S1AB/IMACS_comp_s1ab_ascdesc_iw1_inter.png')
    fig.show()