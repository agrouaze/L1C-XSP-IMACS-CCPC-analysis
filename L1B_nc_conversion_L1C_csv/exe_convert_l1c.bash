#!/bin/bash

/home1/datawork/ljessel/bin/prun --split-max-jobs=300 --mem=500m --name l1cPRUN_S1AB_B07 --max-time=02:00:00 -e /home1/datahome/ljessel/l1c_conversion/convert_l1c.pbs --split-max-lines 100 /home1/datawork/ljessel/l1c_conversion/ALL_S1AB_IW_XSP__1SDV_B07.txt
