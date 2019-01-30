#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("../../python")
from util import hdf5util, windowfunction
import matplotlib.pyplot as plt

PROJECT_DIR = "./"
DATASET_NAME = "apod"

WINDOW_SIZE = 32

#==============================================================================

#apod = windowfunction.get_blackman_window(WINDOW_SIZE)
apod = windowfunction.get_blackman2_window(WINDOW_SIZE)

hdf5util.write_ndarray(apod, PROJECT_DIR + "apod_1d_blackman-" + str(WINDOW_SIZE) + ".h5", DATASET_NAME)

fig = plt.stem(apod)
