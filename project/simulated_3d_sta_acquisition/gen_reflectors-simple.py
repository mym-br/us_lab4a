#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("../../python")

import numpy as np
from util import hdf5util



PROJECT_DIR = "./"
DATASET_NAME = "reflectors"

X = 0.0
X_MIN = -30.0e-3
X_MAX = 30.0e-3
Y = 0.0
Y_MIN = -10.0e-3
Y_MAX = 10.0e-3
Z = 40.0e-3
NUM_POINTS = 7
COEFF_MIN = 0.5
COEFF_MAX = 5.0

def gen_reflectors_file(x, y, z, coeff, n_file):
    reflectors = np.zeros((len(x), 4))
    reflectors[:, 0] = x
    reflectors[:, 1] = y
    reflectors[:, 2] = z
    reflectors[:, 3] = coeff

    hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors" + str(n_file) + ".h5", DATASET_NAME)

gen_reflectors_file(x=[X],
                    y=Y,
                    z=Z,
                    coeff=1.0,
                    n_file=0)

gen_reflectors_file(x=np.linspace(X_MIN, X_MAX, NUM_POINTS),
                    y=Y,
                    z=Z,
                    coeff=1.0,
                    n_file=1)

gen_reflectors_file(x=np.linspace(X_MIN, X_MAX, NUM_POINTS),
                    y=np.linspace(Y_MIN, Y_MAX, NUM_POINTS),
                    z=Z,
                    coeff=1.0,
                    n_file=2)

gen_reflectors_file(x=np.linspace(X_MIN, X_MAX, NUM_POINTS),
                    y=Y,
                    z=Z,
                    coeff=np.linspace(COEFF_MIN, COEFF_MAX, NUM_POINTS),
                    n_file=3)
