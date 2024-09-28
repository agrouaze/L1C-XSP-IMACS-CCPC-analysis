from yaml import CLoader as Loader
import datetime
import pdb
import numpy as np
import os
import l1canalysis
import logging
import glob
import time
from yaml import load
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

sat_colors = {"S1A": "blue", "S1B": "red"}


local_config_pontential_path = os.path.join(
    os.path.dirname(l1canalysis.__file__), "localconfig.yml"
)

if os.path.exists(local_config_pontential_path):
    config_path = local_config_pontential_path
else:
    config_path = os.path.join(os.path.dirname(l1canalysis.__file__), "config.yml")
logging.info("config path: %s", config_path)
stream = open(config_path, "r")
conf = load(stream, Loader=Loader)



def load_and_merge_dataframes(base_path,pola_acqui="1SDV",pola_chosen='VV',burst_type_pattern = 'intraburst',test=False):
    """

    :param base_path: directory where are stored the .csv files example "/home1/scratch/agrouaze/l1c_converted/B07/" # path with only extension
    :param pola: str 1SDV for instance
    :param burst_type_pattern: str intraburst or interburst
    :param pola_chosen: str VV or VH or ...
    :param test: bool True -> reduced number of csv files to load for test/debug
    :return:
    """
    dfs = {}
    for unit in ['S1A','S1B']:
        safe_pattern = unit+'_IW_XSP__'+pola_acqui+'_*.SAFE/'
        logging.info('safe_pattern : %s',safe_pattern)
        #
        # burst_type_pattern = 'interburst'

        ### Only for L1C_4.1 version files
        #csv_polar_pattern = '/*vv*.csv'
        #csv_polar_pattern = '/*vh*.csv'
        #csv_polar_pattern = '/*.csv'

        safes = glob.glob(os.path.join(base_path,safe_pattern,burst_type_pattern))  # we retrieve all the .SAFE files we want
        print('number of SAFE files : '+str(len(safes)))
        fns_csv=[]
        for safe in safes :
            fns_csv += glob.glob(os.path.join(safe,'*.csv'))
            #fns_csv += glob.glob(safe+csv_polar_pattern)      # we retrieve the .csv files from the preselected .SAFE files
        if test is True:
            logging.info('load only 10 .csv files for test')
            fns_csv = fns_csv[0:10]
        print('number of csv files : '+str(len(fns_csv)))


        start_time = time.time()
        dataframes = []                         # Initialise a liste to stock the dataFrames
        for csv_file in fns_csv[:]:             # browse each CSV file
            if os.path.getsize(csv_file) > 0:   # if the file is not empty
                df = pd.read_csv(csv_file)      # read the file and add it to the dataframe
                dataframes.append(df)

        print('number of not empty csv :', len(dataframes))
        df_all = pd.concat(dataframes, ignore_index=True) # Concatenate the dataframes
        end_time = time.time()

        print('loading time :',end_time - start_time,'s')

        df_all = df_all[df_all['land_flag'] == 0]  # we take the tiles that contains only ocean

        ### Selection of the polarisation (Only for L1C_B07 version files)
        # df_vv =
        df_vh = df_all[df_all['pola']=='VH']
        # df = df_vv  # VV polarisation
        dfs['%s_%s'%(unit,'VV')] = df_all[df_all['pola']=='VV']
        dfs['%s_%s'%(unit,'VH')] = df_all[df_all['pola']=='VH']
    df_s1b = dfs['S1B_'+pola_chosen]
    df_s1a = dfs['S1A_'+pola_chosen]
    return df_s1a,df_s1b

def add_ancillary_variables(df_s1a,df_s1b):
    """
    add ancillary wind information

    :param df_s1a:
    :param df_s1b:
    :return:
    """
    dataset = {'S1A':df_s1a,'S1B':df_s1b}
    for sar_unit in dataset:
        u = dataset[sar_unit]['U10']
        v = dataset[sar_unit]['V10']
        dataset[sar_unit]['Wspeed'] = np.sqrt(u**2 + v**2)                                # wind speed norm
        dataset[sar_unit]['wdir'] = wind_direction_from_U_and_V(u_component=u,v_component=v)
        dataset[sar_unit]['wdir_az'] = azimuth_direction(ground_heading_angle=dataset[sar_unit]['ground_heading'],
                                                         u_component=u,v_component=v)
        # df_s1b['wdir'] = np.mod(180+(180/np.pi)*np.arctan2(u,v),360)           # wind direction taken from the north
        # # dataset[sar_unit]['wdir'] = np.mod((180 / np.pi) * np.arctan2(-u, -v), 360) # correction chatgpt
        # dataset[sar_unit]['wdir_az'] = ((dataset[sar_unit]['wdir']-dataset[sar_unit]['ground_heading']-90)%360) # wind direction brought back to the antenna's frame of reference
        # print('-90°')
        # dataset[sar_unit]['wdir_az'] = dataset[sar_unit]['wdir_az'] * -1.
        print('aar version')
        # df_s1b['wdir_az'] = ((df_s1b['wdir'] - df_s1b[
        #     'ground_heading'] ) % 360)  # test agrouaze
        # print('0°') #+90 tested +0 tested -270 tested
        dataset[sar_unit]['sigma0_dB_filt'] = 10*np.log10(dataset[sar_unit]['sigma0_filt'])          # sigma0 filtered in dB


    ### Selection of the variables of interest for the MACS analysis
    # az_wdir_s1b = df_s1b['wdir_az']
    # macs_Im_s1b = df_s1b['macs_Im_lambda_max=50.0']
    # NRCS_s1b = df_s1b['sigma0_dB_filt']
    # incidence_s1b = df_s1b['incidence']
    # wnd_spd_s1b = df_s1b['Wspeed']
    # heading_s1b = df_s1b['ground_heading']
    # descending_s1b = heading_s1b[abs(heading_s1b) > 150]
    # ascending_s1b = heading_s1b[abs(heading_s1b) < 30]

    # ### Computation of others variables from the ones that are in the csv
    # u = df_s1a['U10']
    # v = df_s1a['V10']
    # df_s1a['Wspeed'] = np.sqrt(u**2 + v**2)                                # wind speed norm
    # df_s1a['wdir'] = np.mod(180+(180/np.pi)*np.arctan2(u,v),360)           # wind direction taken from the north
    # df_s1a['wdir_az'] = ((df_s1a['wdir']-df_s1a['ground_heading']-90)%360) # wind direction brought back to the antenna's frame of reference
    # df_s1a['sigma0_dB_filt'] = 10*np.log10(df_s1a['sigma0_filt'])          # sigma0 filtered in dB


    ### Selection of the variables of interest for the MACS/wind analysis
    # az_wdir_s1a = df_s1a['wdir_az']
    # macs_Im_s1a = df_s1a['macs_Im_lambda_max=50.0']
    # NRCS_s1a = df_s1a['sigma0_dB_filt']
    # incidence_s1a = df_s1a['incidence']
    # wnd_spd_s1a = df_s1a['Wspeed']
    # heading_s1a = df_s1a['ground_heading']
    # descending_s1a = heading_s1a[abs(heading_s1a) > 150]
    # ascending_s1a = heading_s1a[abs(heading_s1a) < 30]
    return dataset['S1A'],dataset['S1B']