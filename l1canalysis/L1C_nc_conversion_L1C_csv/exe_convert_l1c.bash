#!/bin/bash

#/home1/datawork/ljessel/bin/prun --split-max-jobs=300 --mem=500m --name l1cPRUN_S1AB_B07 --max-time=02:00:00 -e /home1/datahome/ljessel/l1c_conversion/convert_l1c.pbs --split-max-lines 100 /home1/datawork/ljessel/l1c_conversion/ALL_S1AB_IW_XSP__1SDV_B07.txt
prunexe="/appli/prun/bin/prun"
pbspath="/home1/datahome/agrouaze/git/L1C-XSP_IMACS-analysis/l1canalysis/L1C_nc_conversion_L1C_csv/convert_l1c.pbs"
# voir notebook lister_L1C.ipynb
#listing="/home1/scratch/agrouaze/l1c_conversion/S1AB_IW_L1C_XSP__1SDV_B07_extension_from_creodias_nc_fullpath.txt" # extension dataset
#listing="/home1/scratch/agrouaze/l1c_conversion/S1AB_IW_L1C_XSP__1SDV_B07_medium_plus_extension_nc_fullpath.txt" # dataset entier
#listing="/home1/scratch/agrouaze/l1c_conversion/S1AB_IW_L1C_XSP__1SDV_B07_medium_listing_nc_fullpath.txt" # medium dataset
listing="/home/datawork-cersat-public/project/sarwave/data/listings/L1C_IW_B09_measu_fullpath_ifremer.txt" # training dataset compelte
$prunexe --split-max-jobs=700 --mem=500m --name l1cPRUN_S1AB_B09 --max-time=02:00:00 -e $pbspath $listing
