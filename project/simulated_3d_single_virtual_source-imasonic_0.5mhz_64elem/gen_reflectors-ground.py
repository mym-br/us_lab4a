#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
sys.path.append("../../python")
import numpy as np
from util import hdf5util
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

PROJECT_DIR = "./"
DATASET_NAME = "reflectors"

# m
X_MIN = -0.7
X_MAX =  0.7
Y_MIN = -0.1
Y_MAX =  0.1
Z_MIN =  0.7
Z_MAX =  0.701
NUM_POINTS = 10000
COEFF_MIN = 1.0
COEFF_MAX = 5.0

#==============================================================================

random.seed(42)

reflectors = np.zeros((NUM_POINTS, 4))
for i in range(NUM_POINTS):
    reflectors[i, 0] = random.uniform(X_MIN, X_MAX)
    reflectors[i, 1] = random.uniform(Y_MIN, Y_MAX)
    reflectors[i, 2] = random.uniform(Z_MIN, Z_MAX)
    reflectors[i, 3] = random.uniform(COEFF_MIN, COEFF_MAX)

hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors-ground.h5", DATASET_NAME)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(reflectors[:, 0], reflectors[:, 1], reflectors[:, 2])
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
