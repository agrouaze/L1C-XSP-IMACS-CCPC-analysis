# IMACS analysis by Lucas-Jessel
My codes and notebooks for the anlaysis of Sentinel1 L1B-XSP IMACS parameter made on internship in 2024.

## ***L1B_data_discover*** directory  
- This directory is my begining with L1B data. It contains an example of a sigma0 map (with the Sentinel1 acquisition geometry) over 1 subswath, 3 subswaths and a whole slice. There is also a Notebook to plot a cross-spectra.
- These Notebooks use the map tools stored in the ***map_SAR_variables_toolkit*** directory and taken and adapted from [Lisa Maillard's toolkit](https://gitlab.ifremer.fr/lisa.maillard/my_tools) and also the `one_xspectra_plot.py` script in the ***python_script*** directory to plot the cross-spectra.  

## ***L1C_nc_conversion_L1C_csv*** directory
- Containing all the scripts to convert L1C netcdf (**datasets**) data into L1C csv data (**dataframes**) to make the IMACS analysis easier.
- To convert a large amount of netcdf data into csv I used PRUN, a tool to make "embarasingly parallel" calcul on the IFREMER HPC calculator. The `PRUN_readme.md` explain how to make these massiv conversions step by step. For more information on prun, click [here](https://collab.umr-lops.fr/fr/calcul/calcul/machines/datarmor).


## ***df_generation*** directory 
- Containing 4 Notebooks wich generates dataframes, using the `L1CConverter` in the `/L1C_nc_conversion_L1C_csv/l1c_conversion_approp.py` script, for 1 subswath, 3 subswaths, one day of data and one month of data.
- These Notebooks are a *getting started* of the conversion script `/L1C_nc_conversion_L1C_csv/l1c_conversion_approp.py` before launching the massive conversion on all the data with prun.


## ***python_scripts*** directory 
- containing python scripts required for the IMACS analysis.

## IMACS analysis Notebooks 
- `CMOD_NRCS.ipynb` : Notebook making a plot of NRCS vs Wind direction for a specified wind speed and an incidence using **CMOD5** model
- `Macs_analysis_S1__upgrade_plots.ipynb` : Notebook of the IMACS analysis on S1A and S1B separately. Containing, among other things, the iso-density plots of IMACS vs Wind direction for different wind speed ranges and incidences ranges over the 3 subswaths.
- `MACS_2nd_analysis_S1b.ipynb` : Another IMACS analysis focusing on filtering the data around a GMF of IMACS (not completely finished).
- `S1A_S1B_IMACS_mean_analysis.ipynb` : Notebook of the IMACS analysis wich make an sensor inter-comparison before merging S1B and S1A data to make a complementary analysis of IMACS.
- 


