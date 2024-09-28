from yaml import CLoader as Loader
import datetime
import pdb
import numpy as np
import os
import l1canalysis
import logging
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