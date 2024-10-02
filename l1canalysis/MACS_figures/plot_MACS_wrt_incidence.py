import numpy as np
from matplotlib import pyplot as plt
from l1canalysis.utils import mean_iangle
def macs_wrt_incidence_up_down_cross_wind(df,satellite='S1A+B',part='Re',burstkind='intraburst',polarization='vv',lambda_val='50'):
    ymin = -0.04
    ymax = 0.04
    xmin = 30
    xmax = 46
    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))
    mean_iangle_s1 = list(mean_iangle)
    macs_Im_s1 = df[varname]
    az_wdir = df["wdir_az_scat"]


    wnd_spd_5ms = df["Wspeed"][(df["Wspeed"] > 3) & (df["Wspeed"] < 7)]
    wnd_spd_10ms = df["Wspeed"][(df["Wspeed"] > 8) & (df["Wspeed"] < 12)]
    wnd_spd_15ms = df["Wspeed"][(df["Wspeed"] > 13) & (df["Wspeed"] < 17)]
    wnd_spd_20ms = df["Wspeed"][(df["Wspeed"] > 18) & (df["Wspeed"] < 22)]

    az_up = az_wdir[(az_wdir < 0.5) | (az_wdir > 359.5)];
    az_up_5ms = az_up.loc[az_up.index.intersection(wnd_spd_5ms.index)]
    az_up_10ms = az_up.loc[az_up.index.intersection(wnd_spd_10ms.index)]
    az_up_15ms = az_up.loc[az_up.index.intersection(wnd_spd_15ms.index)]
    az_up_20ms = az_up.loc[az_up.index.intersection(wnd_spd_20ms.index)]
    az_up_tot = [az_up_5ms, az_up_10ms, az_up_15ms, az_up_20ms]

    ### Azimuth = 180° (downwind)
    az_down = az_wdir[(az_wdir <= 180.5) & (az_wdir > 179.5)];
    az_down_5ms = az_down.loc[az_down.index.intersection(wnd_spd_5ms.index)]
    az_down_10ms = az_down.loc[az_down.index.intersection(wnd_spd_10ms.index)]
    az_down_15ms = az_down.loc[az_down.index.intersection(wnd_spd_15ms.index)]
    az_down_20ms = az_down.loc[az_down.index.intersection(wnd_spd_20ms.index)]
    az_down_tot = [az_down_5ms, az_down_10ms, az_down_15ms, az_down_20ms]

    ### Azimuth = 90°  (crosswind)
    az_cross = az_wdir[(abs(az_wdir - 90) < 0.5) | (abs(az_wdir - 270) < 0.5)];
    az_cross_5ms = az_cross.loc[az_cross.index.intersection(wnd_spd_5ms.index)]
    az_cross_10ms = az_cross.loc[az_cross.index.intersection(wnd_spd_10ms.index)]
    az_cross_15ms = az_cross.loc[az_cross.index.intersection(wnd_spd_15ms.index)]
    az_cross_20ms = az_cross.loc[az_cross.index.intersection(wnd_spd_20ms.index)]
    az_cross_tot = [az_cross_5ms, az_cross_10ms, az_cross_15ms, az_cross_20ms]


    fig, ax = plt.subplots(3, 4, figsize=(25, 15))

    ### Azimuth = 0° Upwind
    for i in range(4): # loop over the 4 different configuration of wind speed
        imacs = macs_Im_s1[az_up_tot[i].index]
        Imacs_mean = np.zeros_like(mean_iangle_s1)  # Initialise mean Imacs values
        Imacs_std = np.zeros_like(mean_iangle_s1)  # Initialise mean Imacs values

        for iangle in mean_iangle_s1:
            mask = (df['incidence'][az_up_tot[i].index] >= iangle - 0.4) & (
                        df['incidence'][az_up_tot[i].index] < iangle + 0.4)  # Mask to selects the points in the bin
            ## mean by inc bins
            Imacs_mean[mean_iangle_s1.index(iangle)] = np.mean(imacs[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin
            ## std by inc bins
            Imacs_std[mean_iangle_s1.index(iangle)] = np.std(imacs[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin

        # display(Imacs_mean, Imacs_std)
        ax[0][i].plot(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean)], Imacs_mean[np.isfinite(Imacs_mean)],
                      color='orange', lw=1.5, label='IMACS$_{mean}$', marker='o')
        ax[0][i].fill_between(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean)],
                              Imacs_mean[np.isfinite(Imacs_mean)] - Imacs_std[np.isfinite(Imacs_mean)],
                              Imacs_mean[np.isfinite(Imacs_mean)] + Imacs_std[np.isfinite(Imacs_mean)],
                              label='IMACS$_{STD}$', color='powderblue')

        ax[0][i].hlines(0, 30, 46, color='black', lw=2)
        ax[0][i].tick_params(axis='x', labelsize=10)
        ax[0][i].tick_params(axis='y', labelsize=10)
        ax[0][i].set_xlim(xmin, xmax)
        ax[0][i].set_ylim(ymin, ymax)
        # ax[0][i].set_xlabel('Incidence angle [deg]', fontsize=12)
        if i == 0:
            ax[0][i].set_ylabel(r'$\Im$(MACS) - %s - %sm'%(polarization,lambda_val), fontsize=12)
        else:
            ax[0][i].set_ylabel('')
        ax[0][i].grid(linestyle='--', color='gray', alpha=0.9)
        ax[0][i].set_title(r'wind speed = %d $\pm$ 2 m/s ' % ((i + 1) * 5), fontsize=15)
        lastpart = r'{0;360} $\pm$ 0.5° (upwind)'
        txt_str = 'S1 pts num: %d \nazimuth : \n' % (len(macs_Im_s1[az_up_tot[i].index]))+lastpart
        props = dict(boxstyle='square', facecolor='white')
        ax[0][i].text(0.03, 0.97, txt_str, transform=ax[0][i].transAxes, fontsize=12, verticalalignment='top',
                      bbox=props)
        ax[0][i].legend()

    ### Azimuth = 180°
    # az180 = az_wdir_s1_inter[(az_wdir_s1_inter <= 180.5) & (az_wdir_s1_inter > 179.5)]
    # az180_5ms = az180.loc[az180.index.intersection(wnd_spd_5ms.index)]
    # az180_10ms = az180.loc[az180.index.intersection(wnd_spd_10ms.index)]
    # az180_15ms = az180.loc[az180.index.intersection(wnd_spd_15ms.index)]
    # az180_20ms = az180.loc[az180.index.intersection(wnd_spd_20ms.index)]
    # az180_tot = [az180_5ms, az180_10ms, az180_15ms, az180_20ms]

    for i in range(4):
        imacs = macs_Im_s1[az_down_tot[i].index]
        Imacs_mean = np.zeros_like(mean_iangle_s1)  # Initialise mean Imacs values
        Imacs_std = np.zeros_like(mean_iangle_s1)  # Initialise mean Imacs values

        for iangle in mean_iangle_s1:
            mask = (df['incidence'][az_down_tot[i].index] >= iangle - 0.4) & (
                        df['incidence'][az_down_tot[i].index] < iangle + 0.4)  # Mask to selects the points in the bin
            ## mean by inc bins
            Imacs_mean[mean_iangle_s1.index(iangle)] = np.mean(imacs[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin
            ## std by inc bins
            Imacs_std[mean_iangle_s1.index(iangle)] = np.std(imacs[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin

        # display(Imacs_mean, Imacs_std)
        ax[1][i].plot(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean)], Imacs_mean[np.isfinite(Imacs_mean)],
                      color='orange', lw=1.5, label='IMACS$_{mean}$', marker='o')
        ax[1][i].fill_between(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean)],
                              Imacs_mean[np.isfinite(Imacs_mean)] - Imacs_std[np.isfinite(Imacs_mean)],
                              Imacs_mean[np.isfinite(Imacs_mean)] + Imacs_std[np.isfinite(Imacs_mean)],
                              label='IMACS$_{STD}$', color='powderblue')

        ax[1][i].hlines(0, 30, 46, color='black', lw=2)
        ax[1][i].tick_params(axis='x', labelsize=10)
        ax[1][i].tick_params(axis='y', labelsize=10)
        ax[1][i].set_xlim(xmin, xmax);
        ax[1][i].set_ylim(ymin, ymax)
        if i == 0:
            ax[1][i].set_ylabel(r'$\Im$(MACS) - %s - %sm'%(polarization,lambda_val), fontsize=12)
        else:
            ax[1][i].set_ylabel('')
        ax[1][i].grid(linestyle='--', color='gray', alpha=0.9)
        lastpart = r'%d $\pm$ 0.5° (downwind)'%180
        txt_str = 'S1 pts num: %d \nazimuth : \n' % (
        len(macs_Im_s1[az_down_tot[i].index]))+lastpart
        props = dict(boxstyle='square', facecolor='white')
        ax[1][i].text(0.03, 0.97, txt_str, transform=ax[1][i].transAxes, fontsize=12, verticalalignment='top',
                      bbox=props)
        ax[1][i].legend()

    ### Azimuth = 90°
    # az90 = az_wdir_s1_inter[(abs(az_wdir_s1_inter - 90) < 0.5) | (abs(az_wdir_s1_inter - 270) < 0.5)]
    # az90_5ms = az90.loc[az90.index.intersection(wnd_spd_5ms.index)]
    # az90_10ms = az90.loc[az90.index.intersection(wnd_spd_10ms.index)]
    # az90_15ms = az90.loc[az90.index.intersection(wnd_spd_15ms.index)]
    # az90_20ms = az90.loc[az90.index.intersection(wnd_spd_20ms.index)]
    # az90_tot = [az90_5ms, az90_10ms, az90_15ms, az90_20ms]

    for i in range(4):
        imacs = macs_Im_s1[az_cross_tot[i].index]
        Imacs_mean = np.zeros_like(mean_iangle_s1)  # Initialise mean Imacs values
        Imacs_std = np.zeros_like(mean_iangle_s1)  # Initialise mean Imacs values

        for iangle in mean_iangle_s1:
            mask = (df['incidence'][az_cross_tot[i].index] >= iangle - 0.4) & (
                        df['incidence'][az_cross_tot[i].index] < iangle + 0.4)  # Mask to selects the points in the bin
            ## mean by inc bins
            Imacs_mean[mean_iangle_s1.index(iangle)] = np.mean(imacs[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin
            ## std by inc bins
            Imacs_std[mean_iangle_s1.index(iangle)] = np.std(imacs[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin

        # display(Imacs_mean, Imacs_std)
        ax[2][i].plot(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean)], Imacs_mean[np.isfinite(Imacs_mean)],
                      color='orange', lw=1.5, label='IMACS$_{mean}$', marker='o')
        ax[2][i].fill_between(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean)],
                              Imacs_mean[np.isfinite(Imacs_mean)] - Imacs_std[np.isfinite(Imacs_mean)],
                              Imacs_mean[np.isfinite(Imacs_mean)] + Imacs_std[np.isfinite(Imacs_mean)],
                              label='IMACS$_{STD}$', color='powderblue')

        ax[2][i].hlines(0, 30, 46, color='black', lw=2)
        ax[2][i].tick_params(axis='x', labelsize=10)
        ax[2][i].tick_params(axis='y', labelsize=10)
        ax[2][i].set_xlim(xmin, xmax)
        ax[2][i].set_ylim(ymin, ymax)
        ax[2][i].set_xlabel('Incidence angle [deg]', fontsize=12)
        if i == 0:
            ax[2][i].set_ylabel(r'$\Im$(MACS) - %s - %sm'%(polarization,lambda_val), fontsize=12)
        else:
            ax[2][i].set_ylabel('')
        ax[2][i].grid(linestyle='--', color='gray', alpha=0.9)
        lastpart = r'{90;270} $\pm$ 0.5° (crosswind)'
        txt_str = 'S1 pts num: %d \nazimuth : \n' % (
            len(macs_Im_s1[az_cross_tot[i].index]))+lastpart
        props = dict(boxstyle='square', facecolor='white')
        ax[2][i].text(0.03, 0.97, txt_str, transform=ax[2][i].transAxes, fontsize=12, verticalalignment='top',
                      bbox=props)
        ax[2][i].legend()
    tit = satellite+'IMACS versus incidence angle | wind speed dependency | %s %s'%(burstkind,polarization)
    fig.suptitle(tit, fontsize=15)
    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/s1ab_IMACS_iangle_wnd_filter_az.png')
    fig.show()

def macs_wrt_incidence_recap(df,satellite='S1A+B',part='Re',burstkind='intraburst',polarization='vv',lambda_val='50'):
    ymin = -0.01
    ymax = 0.03

    varname = 'macs_%s_lambda_max=%s' % (part, float(lambda_val))
    mean_iangle_s1 = list(mean_iangle)
    macs_Im_s1 = df[varname]
    az_wdir = df["wdir_az_scat"]

    wnd_spd_5ms = df["Wspeed"][(df["Wspeed"] > 3) & (df["Wspeed"] < 7)]
    wnd_spd_10ms = df["Wspeed"][(df["Wspeed"] > 8) & (df["Wspeed"] < 12)]
    wnd_spd_15ms = df["Wspeed"][(df["Wspeed"] > 13) & (df["Wspeed"] < 17)]
    wnd_spd_20ms = df["Wspeed"][(df["Wspeed"] > 18) & (df["Wspeed"] < 22)]


    # xmin=mean_iangle_s1[0]-0.5;xmax=mean_iangle_s1[12]+0.5
    xmin = mean_iangle_s1[0] - 0.5;
    xmax = mean_iangle_s1[11] + 0.5  # A Grouazel fix
    fig, ax = plt.subplots(1, 4, figsize=(30, 6))

    ### Azimuth = 0/360° (upwind)
    az_up = az_wdir[(az_wdir < 0.5) | (az_wdir > 359.5)]
    az_up_inter = az_up
    az_up_5ms = az_up.loc[az_up.index.intersection(wnd_spd_5ms.index)]
    az_up_10ms = az_up.loc[az_up.index.intersection(wnd_spd_10ms.index)]
    az_up_15ms = az_up.loc[az_up.index.intersection(wnd_spd_15ms.index)]
    az_up_20ms = az_up.loc[az_up.index.intersection(wnd_spd_20ms.index)]
    az_up_tot = [az_up_5ms, az_up_10ms, az_up_15ms, az_up_20ms]

    ### Azimuth = 180° (downwind)
    az_down = az_wdir[(az_wdir <= 180.5) & (az_wdir > 179.5)];
    az_down_inter = az_down
    az_down_5ms = az_down.loc[az_down.index.intersection(wnd_spd_5ms.index)]
    az_down_10ms = az_down.loc[az_down.index.intersection(wnd_spd_10ms.index)]
    az_down_15ms = az_down.loc[az_down.index.intersection(wnd_spd_15ms.index)]
    az_down_20ms = az_down.loc[az_down.index.intersection(wnd_spd_20ms.index)]
    az_down_tot = [az_down_5ms, az_down_10ms, az_down_15ms, az_down_20ms]

    ### Azimuth = 90°  (crosswind)
    az_cross = az_wdir[(abs(az_wdir - 90) < 0.5) | (abs(az_wdir - 270) < 0.5)];
    az_cross_inter = az_cross
    az_cross_5ms = az_cross.loc[az_cross.index.intersection(wnd_spd_5ms.index)]
    az_cross_10ms = az_cross.loc[az_cross.index.intersection(wnd_spd_10ms.index)]
    az_cross_15ms = az_cross.loc[az_cross.index.intersection(wnd_spd_15ms.index)]
    az_cross_20ms = az_cross.loc[az_cross.index.intersection(wnd_spd_20ms.index)]
    az_cross_tot = [az_cross_5ms, az_cross_10ms, az_cross_15ms, az_cross_20ms]

    for i in range(4):
        imacs_cross = macs_Im_s1[az_cross_tot[i].index]
        Imacs_mean_cross = np.zeros_like(mean_iangle_s1)  # Initialise mean crosswind Imacs values
        Imacs_std_cross = np.zeros_like(mean_iangle_s1)  # Initialise mean crosswind Imacs values

        imacs_up = macs_Im_s1[az_up_tot[i].index]
        Imacs_mean_up = np.zeros_like(mean_iangle_s1)  # Initialise mean upwind Imacs values
        Imacs_std_up = np.zeros_like(mean_iangle_s1)  # Initialise mean upwind Imacs values

        imacs_down = macs_Im_s1[az_down_tot[i].index]
        Imacs_mean_down = np.zeros_like(mean_iangle_s1)  # Initialise mean downwind Imacs values
        Imacs_std_down = np.zeros_like(mean_iangle_s1)  # Initialise mean downwind Imacs values

        for iangle in mean_iangle_s1:
            ### upwind
            mask = (df['incidence'][az_up_tot[i].index] >= iangle - 0.4) & (
                        df['incidence'][az_up_tot[i].index] < iangle + 0.4)  # Mask to selects the points in the bin
            ## mean by inc bins
            Imacs_mean_up[mean_iangle_s1.index(iangle)] = np.mean(imacs_up[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin
            ## std by inc bins
            Imacs_std_up[mean_iangle_s1.index(iangle)] = np.std(imacs_up[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs std value in the bin

            ### downwind
            mask = (df['incidence'][az_down_tot[i].index] >= iangle - 0.4) & (df['incidence'][az_down_tot[
                i].index] < iangle + 0.4)  # Mask to selects the points in the bin
            ## mean by inc bins
            Imacs_mean_down[mean_iangle_s1.index(iangle)] = np.mean(imacs_down[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs mean value in the bin
            ## std by inc bins
            Imacs_std_down[mean_iangle_s1.index(iangle)] = np.std(imacs_down[mask]) if np.sum(
                mask) > 0 else np.nan  # Calculation of the Imacs std value in the bin

            ### crosswind
            mask = (df['incidence'][az_cross_tot[i].index] >= iangle - 0.4) & (df['incidence'][az_cross_tot[
                i].index] < iangle + 0.4)  # Mask to selects the points in the bin
            Imacs_mean_cross[mean_iangle_s1.index(iangle)] = np.mean(imacs_cross[mask]) if np.sum(
                mask) > 0 else np.nan  # mean by inc bins
            Imacs_std_cross[mean_iangle_s1.index(iangle)] = np.std(imacs_cross[mask]) if np.sum(
                mask) > 0 else np.nan  # std by inc bins

        ax[i].plot(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean_up)], abs(Imacs_mean_up[np.isfinite(Imacs_mean_up)]),
                   color='blue', lw=1.5, label='upwind', marker='o')
        ax[i].fill_between(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean_up)],
                           abs(Imacs_mean_up[np.isfinite(Imacs_mean_up)]) - Imacs_std_up[np.isfinite(Imacs_mean_up)],
                           abs(Imacs_mean_up[np.isfinite(Imacs_mean_up)]) + Imacs_std_up[np.isfinite(Imacs_mean_up)],
                           color='blue', alpha=0.2)

        ax[i].plot(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean_down)],
                   abs(Imacs_mean_down[np.isfinite(Imacs_mean_down)]), color='red', lw=1.5, label='downwind',
                   marker='o')
        ax[i].fill_between(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean_down)],
                           abs(Imacs_mean_down[np.isfinite(Imacs_mean_down)]) - Imacs_std_down[
                               np.isfinite(Imacs_mean_down)],
                           abs(Imacs_mean_down[np.isfinite(Imacs_mean_down)]) + Imacs_std_down[
                               np.isfinite(Imacs_mean_down)],
                           color='red', alpha=0.2)

        ax[i].plot(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean_cross)],
                   abs(Imacs_mean_cross[np.isfinite(Imacs_mean_cross)]), color='green', lw=1.5, label='crosswind',
                   marker='o')
        ax[i].fill_between(np.array(mean_iangle_s1)[np.isfinite(Imacs_mean_cross)],
                           abs(Imacs_mean_cross[np.isfinite(Imacs_mean_cross)]) - Imacs_std_cross[
                               np.isfinite(Imacs_mean_cross)],
                           abs(Imacs_mean_cross[np.isfinite(Imacs_mean_cross)]) + Imacs_std_cross[
                               np.isfinite(Imacs_mean_cross)],
                           color='green', alpha=0.2)

        ax[i].hlines(0, 30, 46, color='black', lw=2)
        ax[i].tick_params(axis='x', labelsize=15)
        ax[i].tick_params(axis='y', labelsize=15)
        ax[i].set_xlim(xmin, xmax)
        ax[i].set_ylim(ymin, ymax)
        ax[i].set_xlabel('Incidence angle [deg]', fontsize=15)
        if i == 0:
            ax[i].set_ylabel('|$\Im$(MACS)| - %s - %sm'%(polarization,lambda_val), fontsize=15)
        else:
            ax[i].set_ylabel('')
        ax[i].grid(linestyle='--', color='gray', alpha=0.9)
        ax[i].set_title('wind speed = %d $\pm$ 2 m/s ' % ((i + 1) * 5), fontsize=15)
        txt_str = 'S1 pts num: %d' % (
                    len(macs_Im_s1[az_up_tot[i].index]) + len(macs_Im_s1[az_cross_tot[i].index]) + len(
                macs_Im_s1[az_down_tot[i].index]))
        props = dict(boxstyle='square', facecolor='white')
        ax[i].text(0.03, 0.97, txt_str, transform=ax[i].transAxes, fontsize=15, verticalalignment='top', bbox=props)
        ax[i].legend(loc='lower center', fontsize=10)

    fig.suptitle(satellite+' | IMACS versus incidence angle | wind speed dependency | %s %s' % (burstkind,polarization), fontsize=15)
    # fig.savefig('/home1/datahome/ljessel/Plots/MACS_analysis/inter/S1A+B/IMACSvsiangle_wnd_filter_s1ab_allIW_inter.png')
    fig.show()