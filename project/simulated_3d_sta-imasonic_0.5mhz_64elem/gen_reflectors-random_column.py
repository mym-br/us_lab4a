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

X_CENTER = 0.2
Y_CENTER = 0.0
Z_MIN = 0.6
Z_MAX = 1.0
MAX_RADIUS = 30.0e-3
NUM_POINTS = 30
COEFF_MIN = 1.0
COEFF_MAX = 5.0

#==============================================================================

random.seed(42)

reflectors = np.zeros((NUM_POINTS, 4))
for i in range(NUM_POINTS):
    r = abs(random.gauss(0.0, MAX_RADIUS / 3.0))
    theta = random.uniform(0.0, 2.0 * np.pi)
    z0 = (Z_MIN + Z_MAX) / 2.0
    z = abs(random.gauss(z0, 0.5 * (Z_MAX - Z_MIN) / 3.0))
    reflectors[i, 0] = X_CENTER + r * np.cos(theta)
    reflectors[i, 1] = Y_CENTER + r * np.sin(theta)
    reflectors[i, 2] = z
    reflectors[i, 3] = COEFF_MIN + (COEFF_MAX - COEFF_MIN) * random.random()

hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors-random_column.h5", DATASET_NAME)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(reflectors[:, 0], reflectors[:, 1], reflectors[:, 2])
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
