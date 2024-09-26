"""ASA_WVI_1PNP
Developped for the version 4.1 of l1c processor.
"""

import os
import glob
import pandas as pd
import xarray as xr
import logging

class L1CConverter:
    def __init__(self, path_l1c_subswath, root_savepath=None, selected_vars=['sigma0', 'incidence', 'macs_Im','macs_Re'], burst_type='intraburst',polarisations=['VV','VH']):
        """
        Initialize the L1C_Converter.

        Parameters:
        - path_l1c (str): Path to the l1c subswath file.
        - root_savepath (str): Root path where the converted file will be saved.
        - selected_vars (list): List of selected variables to extract from the dataset.
        - burst_type (str): Either 'interburst' or 'intraburst'.
        """
        self.path = path_l1c_subswath
        self.root_savepath = root_savepath
        self.selected_vars = selected_vars
        self.burst_type = burst_type
        self.polarisations = polarisations

        
    def converter(self, save=True):
        """
        Convert the l1c subswath file to a dataframe.
        
        Parameters
        - save (boolean): If True the converted file is saved to a CSV. If False, the converted dataframe is returned. Defaults to True.

        Returns:
        - None if the file cannot be converted. Else, see save parameter.
        """

        ds = xr.open_dataset(self.path, group=self.burst_type)
        res = self.get_dataframe(ds)

        if res is not None:
            
            if save:
                output = self.save_dataframe(res)
                return output
                
            else:
                return res

    # def get_df_cwaves(self,ds):
    #     logging.debug('ds : %s',ds['cwave_params'])
    #     df_cwaves = pd.DataFrame()
    #     for pol in self.polarisations:
    #         for _phi_hf in ds['phi_hf'].values:
    #             for _k_gp in ds['k_gp'].values:
    #                 varname = f'cwave_params_pol={pol}_k_gp={_k_gp}_and_phi_hf={_phi_hf}'
    #                 df_cwaves[varname] = ds['cwave_params'].sel(k_gp=_k_gp,
    #                     phi_hf=_phi_hf).to_dataframe().reset_index(drop=True).drop(columns=['k_gp', 'phi_hf', 'spatial_ref', 'sample', 'line', 'pol'], errors='ignore')
    #     return df_cwaves

    def get_df_imacs(self,ds):
        df_imacs = None
        for lambdav in ds['lambda_range_max_macs'].values:
            for pol in self.polarisations:
                varname = f'imacs_re_l={lambdav}_pol={pol}'
                tmp = ds['macs_Re'].sel(pol=pol,
                    lambda_range_max_macs=lambdav).to_dataframe().reset_index(drop=True).drop(columns=['spatial_ref', 'sample', 'line','pol','lambda_range_max_macs'], errors='ignore')
                tmp.rename(columns={'macs_Re':varname},inplace=True,errors='raise')
                # logging.debug('tmp : %s',tmp)
                if df_imacs is None:
                    df_imacs = tmp
                else:
                    df_imacs = df_imacs.merge(tmp,on=['longitude', 'latitude'])
                # df_imacs[varname] = tmp
                varname = f'imacs_im_l={lambdav}_pol={pol}'
                tmp = ds['macs_Im'].sel(pol=pol,
                    lambda_range_max_macs=lambdav).to_dataframe().reset_index(drop=True).drop(columns=['spatial_ref', 'sample', 'line','pol','lambda_range_max_macs'], errors='ignore')
                tmp.rename(columns={'macs_Im':varname},inplace=True,errors='raise')
                df_imacs = df_imacs.merge(tmp,on=['longitude', 'latitude'])
        return df_imacs
    
    def get_dataframe(self, ds):
        """
        Convert Sentinel-1 Level-1C sub-swath data to a pandas DataFrame.

        Parameters:
        - ds (xarray.Dataset): Input Sentinel-1 Level-1C sub-swath dataset.
        - path (str): File path corresponding to the dataset.

        Returns:
        - pd.DataFrame : Returns a pandas DataFrame containing the extracted information from the input dataset.
        """
        if not set(self.selected_vars).issubset(ds.keys()):
            print(f'All Variables not found in {self.path}')
            for vv in self.selected_vars:
                if vv not in ds:
                    logging.debug('absent variable: %s',vv)
            return None
            
        ds = ds.drop_vars('crs')
        ds = ds[self.selected_vars]
        logging.debug('ds for selected vars: %s',ds)
        # a base dataframe with only variables without polarization
        # df_res = ds.drop_dims(['k_gp', 'phi_hf','pol']).squeeze().to_dataframe().reset_index(drop=True).drop(columns=['spatial_ref',
        #     'tile_sample', 'tile_line'], errors='ignore')
        basic_variables = [vv for vv in ds if 'pol' not in ds[vv].coords]
        df_res = ds[basic_variables].squeeze().to_dataframe().reset_index(drop=True).drop(columns=['spatial_ref','sample', 'line'], errors='ignore')
        # df_polarized = pd.DataFrame()
        df_polarized = None

        for vv in ds.variables:
            logging.debug('\n ========== var %s %s',vv,type(vv))
            if 'pol' in ds[vv].coords and vv not in ['cwave_params','macs_Im','macs_Re','pol','spatial_ref']:
                for polval in self.polarisations:
                    tmp = ds[vv].sel(pol=polval).to_dataframe().reset_index(drop=True).drop(columns=['spatial_ref','sample', 'line','pol'], errors='ignore')
                    # logging.debug('tmp : %s',tmp)
                    varnamee = vv+'_'+polval
                    tmp.rename(columns={vv:varnamee},inplace=True,errors='raise')
                    if df_polarized is None:
                        # logging.debug('first definition of df_polarized')
                        df_polarized = tmp
                    else:
                        # logging.debug('merge of df_polarized')
                        df_polarized = df_polarized.merge(tmp,on=['longitude', 'latitude'])
                    # logging.debug('%s',df_polarized.keys())
                    # logging.debug('polarized %s %s',vv+'_'+polval,df_polarized[varnamee])
        # df_cwaves = self.get_df_cwaves(ds)
        df_imacs = self.get_df_imacs(ds)
        df_res = df_res.merge(df_imacs, on=['longitude', 'latitude'])
        # df_res = df_res.merge(df_cwaves, on=['longitude', 'latitude'])
        df_res = df_res.merge(df_polarized,on=['longitude', 'latitude'])


        df_res = df_res.assign(file_path=self.path)
        return df_res
    
        
    def save_dataframe(self, df):
        """
        Save the converted data to a CSV file.

        Parameters:
        - df (pandas.DataFrame): Converted l1c subswath data.

        The CSV file is saved with the same base name as the corresponding NetCDF file.
        """
        
        safe_name, file_name = self.path.split(os.sep)[-2:]
        savepath = os.path.join(self.root_savepath, safe_name, self.burst_type, file_name).replace('.nc', '.csv')

        os.makedirs(os.path.dirname(savepath), exist_ok=True)
        df.to_csv(savepath, index=False)
        return savepath