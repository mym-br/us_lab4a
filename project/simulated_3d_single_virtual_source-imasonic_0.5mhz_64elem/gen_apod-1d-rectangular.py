#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("../../python")
import numpy as np
from util import hdf5util
import matplotlib.pyplot as plt

PROJECT_DIR = "./"
DATASET_NAME = "apod"

WINDOW_SIZE = 32

#==============================================================================

apod = np.ones(WINDOW_SIZE)

hdf5util.write_ndarray(apod, PROJECT_DIR + "apod_1d_rectangular-" + str(WINDOW_SIZE) + ".h5", DATASET_NAME)

fig = plt.stem(apod)
