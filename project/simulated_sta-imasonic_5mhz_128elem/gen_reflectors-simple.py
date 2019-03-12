#!/usr/bin/env python3
# This file is in the public domain.

import sys
sys.path.append("../../python")

import numpy as np
from util import hdf5util



PROJECT_DIR = "./"
DATASET_NAME = "reflectors"

X = 0.0
X_MIN = -15.0e-3
X_MAX = 15.0e-3
Z = 20.0e-3
NUM_POINTS = 5

def gen_reflectors_file(x, z, n_file):
    reflectors = np.zeros((len(x), 2))
    reflectors[:, 0] = x
    reflectors[:, 1] = z

    hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors" + str(n_file) + ".h5", DATASET_NAME)

gen_reflectors_file(x=[X],
                    z=Z,
                    n_file=0)

gen_reflectors_file(x=np.linspace(X_MIN, X_MAX, NUM_POINTS),
                    z=Z,
                    n_file=1)
