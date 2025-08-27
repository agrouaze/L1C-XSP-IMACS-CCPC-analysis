from matplotlib import pyplot as plt
import numpy as np
import logging
import os
import pickle
import pandas as pd
import l1canalysis
from l1canalysis.MACS_figures import plot_MACS_wrt_incidence
from l1canalysis.MACS_figures import plot_MACS_wrt_windspeed
from l1canalysis.MACS_figures import plot_MACS_azimuth_modulation

def imacs_vs_azimuth_per_windspeed(imacs_interp,imacsformula,fout):
    phi = np.deg2rad(np.arange(0, 360, 1))
    plt.figure(figsize=(16, 12), dpi=100)
    for wsi, u in enumerate([6, 10, 15, 19]):
        plt.subplot(2, 2, wsi + 1)
        plt.plot(np.arange(0, 360, 1), imacsformula(phi, 34, u), 'g', label='IMACS prediction 34°')
        plt.plot(np.arange(0, 360, 1), imacsformula(phi, 42, u), ':g', label='IMACS prediction 42°')
        # plt.plot(np.arange(0, 360, 1), imacs_fit15_agrouaze(phi, 34, u), 'b', label='NN fit agrouaze15 34°')
        # plt.plot(np.arange(0, 360, 1), imacs_fit15_agrouaze(phi, 42, u), ':b', label='NN fit agrouaze15 42°')
        res_org1 = imacs_interp((u, np.arange(0, 360, 1), 34))
        res_org2 = imacs_interp((u, np.arange(0, 360, 1), 42))
        plt.plot(np.arange(0, 360, 1), res_org1, 'r', label='interpolator 34°', lw=3, alpha=0.5)
        plt.plot(np.arange(0, 360, 1), res_org2, ':r', label='interpolator 42°', lw=3, alpha=0.5)
        plt.legend(bbox_to_anchor=(1, 1))
        plt.xlabel("Azimuth, deg")
        plt.ylabel("IMACS")
        plt.title("Wind speed = " + str(u) + 'm/s')
        plt.grid(True)

    plt.savefig(fout)
    logging.info(' output figure: %s',fout)
    # plt.show()

def imacs_vs_incidence(imacs_interp,imacsformula,fout):

    # plt.close('all')
    plt.figure(figsize=(8, 6), dpi=100)
    phi = np.deg2rad(0)
    alpha = np.arange(32, 46, 1)
    alr = np.deg2rad(alpha)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    for u in np.arange(6, 20, 2):
        ff3 = imacsformula(alpha=alpha, u=u, phi=phi)
        rr = imacs_interp((u, np.rad2deg(phi), alpha))
        # choose a random color for the plot color
        # color = random.choice(colors)
        # plt.plot(alpha, ff, '-', color=color, label="fit 1 " + str(u))
        line1, = plt.plot(alpha, ff3, '-', label="analytical " + str(u) + ' m/s')
        line1_color = line1.get_color()
        plt.plot(alpha, rr, ":", color=line1_color, label="interpolator " + str(u) + ' m/s')
    plt.legend(ncols=2, bbox_to_anchor=(1, 1))
    plt.xlabel("Incidence angle, °")
    plt.ylabel("IMACS []")
    plt.title("phi = %i°" % np.rad2deg(phi))
    plt.grid(True)
    plt.savefig(fout)
    logging.info(' output figure: %s', fout)

def imacs_vs_inc_12subplots(imacs_interp,my_pola,lambda_val,burst_family,fout):
    product_id = 'B09'  # it is basically the SARWAVE training dataset
    data_dir = '/home/datawork-cersat-public/cache/project/sarwave/data/products/developments/L1C_dataframes/slc/iw/'
    f_temporary = os.path.join(data_dir, product_id,
                               'temporary_df_iw_l1c_macs_ccpc_%s_product_%s_test_calib_windspeed0windir0_%s.pkl' % (
                               burst_family, product_id, my_pola))
    fid = open(f_temporary, 'rb')
    tmpdata = pickle.load(fid)
    fid.close()

    maxvava = 0.9
    fig, ax = l1canalysis.MACS_figures.plot_MACS_wrt_incidence.macs_wrt_incidence_up_down_cross_wind(
        tmpdata['dfs_l1c']['S1A_' + my_pola],
        satellite='S1A+B', part='Im',
        burstkind=burst_family, polarization=my_pola.lower(), lambda_val=lambda_val, ymax=maxvava)
    # ax[0,0].plot(38,0,'ro',ms=8)
    # add fit curves
    # fig, ax = plt.subplots(3, 4, figsize=(25, 15))
    azicases = {'upwind': 0, 'downwind': 180, 'crosswind': 90}
    incidence_vect = np.arange(32, 45, 0.5)
    for plotx, aziX in enumerate(azicases):
        for ploty, windK in enumerate([5, 10, 15, 20]):
            # print(azicases[aziX],windK)
            try:
                rr = imacs_interp((windK, azicases[aziX], incidence_vect))
                ax[plotx, ploty].plot(incidence_vect, rr, '.-', lw=3, ms=4, markeredgecolor='k',
                                      label='interpolator')
            except:
                print('error', azicases[aziX], windK)
                pass
            ax[plotx, ploty].set_ylim(-maxvava, maxvava)
            ax[plotx, ploty].grid(True)
            ax[plotx, ploty].axhline(y=0, c='k', alpha=0.5)
            ax[plotx, ploty].tick_params(axis='x', labelsize=10)
            ax[plotx, ploty].tick_params(axis='y', labelsize=10)
            if ploty == 0:
                ax[plotx, ploty].set_ylabel(r'$\Im$(MACS) - %s - %sm' % (my_pola, lambda_val), fontsize=12)
            else:
                ax[plotx, ploty].set_ylabel('')
            ax[plotx, ploty].grid(linestyle='--', color='gray', alpha=0.9)
            ax[plotx, ploty].set_title(r'wind speed = %d $\pm$ 2 m/s ' % (windK), fontsize=15)
            if plotx == 2:
                ax[plotx, ploty].set_xlabel('Incidence angle [deg]', fontsize=12)
    # plt.show(fig)
    plt.savefig(fout)
    logging.info(' output figure: %s', fout)

def imacs_vs_windspeed_fourches(imacsformula,my_pola,lambda_val,burst_family,sar_unit,fout):
    azicases = {'upwind': 0, 'downwind': 180, 'crosswind': 90}
    product_id = 'B09'  # it is basically the SARWAVE training dataset
    data_dir = '/home/datawork-cersat-public/cache/project/sarwave/data/products/developments/L1C_dataframes/slc/iw/'
    f_temporary = os.path.join(data_dir, product_id,
                               'temporary_df_iw_l1c_macs_ccpc_%s_product_%s_test_calib_windspeed0windir0_%s.pkl' % (
                                   burst_family, product_id, my_pola))
    fid = open(f_temporary, 'rb')
    tmpdata = pickle.load(fid)
    fid.close()

    maxvava = 0.9
    if True:
        fig, ax, incs = l1canalysis.MACS_figures.plot_MACS_wrt_windspeed.fig_MACS_row_incidence_multi_azi_wrt_windspeed(
            tmpdata['dfs_l1c'][sar_unit + '_' + my_pola],
            satellite=sar_unit, part='Im',
            burstkind=burst_family,
            polarization=my_pola,
            lambda_val=lambda_val, ymax=maxvava)

    windspeedvect = np.arange(3, 21, 0.1)
    windspeedvect_large = np.arange(0, 30, 0.5)
    colourrz = ['#1f77b4', '#ff7f0e', '#2ca02c']  # default bleu orange vert
    # fig, ax = plt.subplots(len(incs),1, figsize=(15, 15))
    for ii, inc in enumerate(incs):
        for plotx, aziX in enumerate(azicases):
            # rr = imacs_interp((windspeedvect, azicases[aziX], inc))
            # ax[ii].plot(windspeedvect,rr,'.-',lw=3,ms=4,markeredgecolor='k',label='interpolator %s'%aziX)
            af = imacsformula(alpha=inc, u=windspeedvect_large, phi=np.deg2rad(azicases[aziX]))
            ax[ii].plot(windspeedvect_large, af, '--', c=colourrz[plotx], lw=3, ms=4, markeredgecolor='k',
                        label='analytical form %s' % aziX)
        ax[ii].legend(bbox_to_anchor=(1, 1))
    plt.savefig(fout)
    logging.info(' output figure: %s', fout)


def imacs_vs_azimuth_4subplots(imacs_interp,imacsformula,my_pola,lambda_val,burst_family,sar_unit,fout):
    maxvava = 0.9
    product_id = 'B09'  # it is basically the SARWAVE training dataset
    data_dir = '/home/datawork-cersat-public/cache/project/sarwave/data/products/developments/L1C_dataframes/slc/iw/'
    f_temporary = os.path.join(data_dir, product_id,
                               'temporary_df_iw_l1c_macs_ccpc_%s_product_%s_test_calib_windspeed0windir0_%s.pkl' % (
                                   burst_family, product_id, my_pola))
    fid = open(f_temporary, 'rb')
    tmpdata = pickle.load(fid)
    fid.close()

    fig, ax, vect_incidence, colors, obs_imacs_mean = l1canalysis.MACS_figures.plot_MACS_azimuth_modulation.macs_az_windspeed_inc_recap(
        tmpdata['dfs_l1c'][sar_unit + '_' + my_pola],
        satellite=sar_unit, part='Im', burstkind=burst_family,
        polarization=my_pola, lambda_val=lambda_val, ymax=maxvava)

    # vect_azimuth = np.arange(0,360,1)
    delta_azi = 5
    vect_azimuth = np.arange(delta_azi / 2, 360 + delta_azi / 2,
                             delta_azi)  # same definition than in utils.py|mean_curve_calc2()

    # vect_azimuth = np.array([1.4,0.2])
    print(vect_incidence)
    # vect_incidence = np.array([34,41])
    # vect_incidence = vect_incidence[1:-1]
    # colors = colors[1:-1]
    # vect_incidence_filtered = vect_incidence[vect_incidence != 36.79411497]
    interpolator_imacs_mean = {}
    analyticalform_imacs_mean = {}
    displayed_var = 'analyticalform'
    for ploty, windK in enumerate([5, 10, 15, 20]):
        for plotx, inc in enumerate(vect_incidence):
            rr = imacs_interp((windK, vect_azimuth, inc))
            interpolator_imacs_mean['ws%s_inc%s' % (windK, inc)] = rr
            # profile_imacs_af = imacs_fit14_agrouaze(alpha=inc, u=windK, phi=np.deg2rad(vect_azimuth))
            profile_imacs_af = imacsformula(alpha=inc, u=windK, phi=np.deg2rad(vect_azimuth))
            analyticalform_imacs_mean['ws%s_inc%s' % (windK, inc)] = profile_imacs_af
            if displayed_var == 'analyticalform':
                ax[ploty].plot(vect_azimuth, profile_imacs_af, '--', lw=3, ms=4, markeredgecolor='k',
                               label='analytical form inc: %1.1f°' % inc, color=colors[plotx])
            elif displayed_var == 'interpolator':
                ax[ploty].plot(vect_azimuth, rr, '--', lw=3, ms=4, markeredgecolor='k',
                               label='interpolator inc: %1.1f°' % inc, color=colors[plotx])
            else:
                raise ValueError('%s'%displayed_var)
        ax[ploty].legend()
    plt.savefig(fout)
    logging.info(' output figure: %s', fout)
    return analyticalform_imacs_mean,interpolator_imacs_mean,obs_imacs_mean

def get_table_residuals_obs_vs_interpoaltor(obs_imacs_mean,interpolator_imacs_mean):
    obs_imacs_mean.keys()
    residus = pd.DataFrame()
    residus_mean = []
    residus_std = []
    mean_obs = []
    mean_interp = []
    index = []
    for kk in obs_imacs_mean:
        index.append(kk)
        xxx = obs_imacs_mean[kk]
        yyy = interpolator_imacs_mean[kk]
        diffz = xxx - yyy
        residus_mean.append(np.mean(diffz))
        residus_std.append(np.std(diffz))
        mean_obs.append(np.mean(xxx))
        mean_interp.append(np.mean(yyy))
    residus['residus_mean'] = residus_mean
    residus['residus_std'] = residus_std
    residus['mean_obs'] = mean_obs
    residus['mean_interp'] = mean_interp
    residus.index = index
    return residus


def get_table_residuals_obs_vs_analyticalform(obs_imacs_mean,analyticalform_imacs_mean):
    obs_imacs_mean.keys()
    residus = pd.DataFrame()
    residus_mean = []
    residus_std = []
    mean_obs = []
    mean_interp = []
    index = []
    for kk in obs_imacs_mean:
        index.append(kk)
        xxx = obs_imacs_mean[kk]
        # yyy = interpolator_imacs_mean[kk]
        yyy = analyticalform_imacs_mean[kk]
        diffz = xxx - yyy
        residus_mean.append(np.mean(diffz))
        residus_std.append(np.std(diffz))
        mean_obs.append(np.mean(xxx))
        mean_interp.append(np.mean(yyy))
    residus['residus_mean'] = residus_mean
    residus['residus_std'] = residus_std
    residus['mean_obs'] = mean_obs
    residus['mean_analytical_form'] = mean_interp
    residus.index = index
    return residus