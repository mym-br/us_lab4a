#!/usr/bin/env python3
# This file is in the public domain.

import sys
sys.path.append("../../../python")
import numpy as np
from util import hdf5util
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

PROJECT_DIR = "./"
DATASET_NAME = "reflectors"

# m
GROUND_X_MIN = -0.7
GROUND_X_MAX =  0.7
GROUND_Y_MIN = -0.1
GROUND_Y_MAX =  0.1
GROUND_Z_MIN =  0.7
GROUND_Z_MAX =  0.701

GROUND_NUM_POINTS = 1000
GROUND_COEFF_MIN = 1.0
GROUND_COEFF_MAX = 5.0

NUM_POINT_PAIRS = 3
# m
PAIR_DISTANCE = 0.05
PAIR_X_MIN = -0.1
PAIR_X_MAX = 0.1
PAIR_Y_MIN = 0.0
PAIR_Y_MAX = 0.0
PAIR_Z_MIN = 0.3
PAIR_Z_MAX = 0.5
PAIR_COEFF = 1.0

#==============================================================================

random.seed(42)

reflectors = np.zeros((GROUND_NUM_POINTS + 2 * NUM_POINT_PAIRS, 4))

for i in range(GROUND_NUM_POINTS):
    reflectors[i, 0] = random.uniform(GROUND_X_MIN, GROUND_X_MAX)
    reflectors[i, 1] = random.uniform(GROUND_Y_MIN, GROUND_Y_MAX)
    reflectors[i, 2] = random.uniform(GROUND_Z_MIN, GROUND_Z_MAX)
    reflectors[i, 3] = random.uniform(GROUND_COEFF_MIN, GROUND_COEFF_MAX)

dx = (PAIR_X_MAX - PAIR_X_MIN) / (NUM_POINT_PAIRS - 1)
dy = (PAIR_Y_MAX - PAIR_Y_MIN) / (NUM_POINT_PAIRS - 1)
dz = (PAIR_Z_MAX - PAIR_Z_MIN) / (NUM_POINT_PAIRS - 1)
for i in range(NUM_POINT_PAIRS):
    offset = GROUND_NUM_POINTS
    reflectors[offset + i, 0] = PAIR_X_MIN + i * dx - 0.5 * PAIR_DISTANCE
    reflectors[offset + i, 1] = PAIR_Y_MIN + i * dy
    reflectors[offset + i, 2] = PAIR_Z_MIN + i * dz
    reflectors[offset + i, 3] = PAIR_COEFF
for i in range(NUM_POINT_PAIRS):
    offset = GROUND_NUM_POINTS + NUM_POINT_PAIRS
    reflectors[offset + i, 0] = PAIR_X_MIN + i * dx + 0.5 * PAIR_DISTANCE
    reflectors[offset + i, 1] = PAIR_Y_MIN + i * dy
    reflectors[offset + i, 2] = PAIR_Z_MIN + i * dz
    reflectors[offset + i, 3] = PAIR_COEFF

hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors-point_pairs_and_ground.h5", DATASET_NAME)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(reflectors[:, 0], reflectors[:, 1], reflectors[:, 2])
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
