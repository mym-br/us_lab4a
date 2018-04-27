#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from util import hdf5util



PROJECT_DIR = "../lab4a/simulated_3d_sta_acquisition/"
DATASET_NAME = "reflectors"

def gen_reflectors_file(x, y, z, coeff, n_file):
    reflectors = np.zeros((len(x), 4))
    reflectors[:, 0] = x
    reflectors[:, 1] = y
    reflectors[:, 2] = z
    reflectors[:, 3] = coeff

    hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors" + str(n_file) + ".h5", DATASET_NAME)

gen_reflectors_file(x=[0.0],
                    y=0.0,
                    z=40.0e-3,
                    coeff=1.0,
                    n_file=0)

gen_reflectors_file(x=np.linspace(-30.0e-3, 30.0e-3, 7),
                    y=0.0,
                    z=40.0e-3,
                    coeff=1.0,
                    n_file=1)

gen_reflectors_file(x=np.linspace(-30.0e-3, 30.0e-3, 7),
                    y=np.linspace(-10.0e-3, 10.0e-3, 7),
                    z=40.0e-3,
                    coeff=1.0,
                    n_file=2)

gen_reflectors_file(x=np.linspace(-30.0e-3, 30.0e-3, 7),
                    y=0.0,
                    z=40.0e-3,
                    coeff=np.linspace(0.5, 5.0, 7),
                    n_file=3)
