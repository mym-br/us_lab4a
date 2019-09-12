#!/usr/bin/env python3
# This file is in the public domain.

import sys
sys.path.append("../../../python")

import numpy as np
from util import hdf5util



PROJECT_DIR = "./"
DATASET_NAME = "reflectors"

# x = z * z_factor
X_Z_FACTORS = [0.0, 0.3, 0.6]

Y = 0.0
Z_MIN = 100.0e-3
Z_MAX = 900.0e-3
Z_LEVELS = 9
COEFF = 1.0

def gen_reflectors_file(x, y, z, coeff):
    reflectors = np.zeros((len(x), 4))
    reflectors[:, 0] = x
    reflectors[:, 1] = y
    reflectors[:, 2] = z
    reflectors[:, 3] = coeff

    hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors-fan.h5", DATASET_NAME)

x_list = []
z_list = []
for point_z in np.linspace(Z_MIN, Z_MAX, Z_LEVELS):
    for z_factor in X_Z_FACTORS:
        x_list.append(point_z * z_factor)
        z_list.append(point_z)

gen_reflectors_file(x=x_list,
                    y=Y,
                    z=z_list,
                    coeff=COEFF)
