"""ASA_WVI_1PNP
Developped for the version B07 of l1c processor.
"""

import os
import glob
import pandas as pd
import xarray as xr
import numpy as np


class L1CConverter:
    def __init__(self, path_l1c_subswath, root_savepath=None, selected_vars=['sigma0', 'incidence', 'macs_Im','macs_Re'], burst_type='intraburst'):
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
                self.save_dataframe(res)
                
            else:
                return res

    
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
            return None
            
        ds = ds.drop_vars('crs')
        ds = ds[self.selected_vars]
                
        ### If we consider the lambda_max as a vector, we need to iterate
        
        ### For the processor B07 : ###
        tile_line = ds.sizes['tile_line']; tile_sample = ds.sizes['tile_sample']
        pol_values = ds['pol'].values  # ['VH', 'VV']
        tile_line_indices, tile_sample_indices, pol_indices = np.meshgrid(np.arange(tile_line), np.arange(tile_sample), np.arange(len(pol_values)), indexing='ij')  # make a grid for the polarisation dimension
        pol_indices = pol_indices.ravel()   # Reorganise the indices so that they are in linear form        
        expanded_pol = pol_values[pol_indices]  # obtain the 'pol' values associated to the indices
        
        df_res = ds.drop_dims(['lambda_range_max_macs']).squeeze().to_dataframe().reset_index(drop=True).drop(columns=['spatial_ref', 'sample', 'line', 'pol'], errors='ignore') # convert the dataset in a dataframe
        df_res['pola'] = expanded_pol      # add 'pola' column with the good values associated

        for _lambda_max in ds['lambda_range_max_macs'].values:
            df = ds[['macs_Im', 'macs_Re']].sel(lambda_range_max_macs = _lambda_max).to_dataframe().reset_index(drop=True).drop(columns=['lambda_range_max_macs', 'spatial_ref', 'sample', 'line', 'pol'], errors='ignore')
            df = df.rename(columns={'macs_Im': f'macs_Im_lambda_max={_lambda_max}',
                                    'macs_Re': f'macs_Re_lambda_max={_lambda_max}'})
            df['pola']=expanded_pol
            df_res = df_res.merge(df, on=['longitude', 'latitude','pola'])


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
