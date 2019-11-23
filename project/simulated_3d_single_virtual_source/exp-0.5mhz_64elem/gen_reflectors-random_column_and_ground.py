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
GROUND_X_MIN = -0.4
GROUND_X_MAX =  0.4
GROUND_Y_MIN = -0.1
GROUND_Y_MAX =  0.1
GROUND_Z_MIN =  0.35
GROUND_Z_MAX =  0.351

GROUND_NUM_POINTS = 1000
GROUND_COEFF_MIN = 1.0
GROUND_COEFF_MAX = 5.0

COLUMN_NUM_POINTS = 30
COLUMN_X_CENTER = 0.1
COLUMN_Y_CENTER = 0.0
COLUMN_Z_MIN = 0.15
COLUMN_Z_MAX = 0.35
COLUMN_MAX_RADIUS = 15.0e-3
COLUMN_COEFF_MIN = 1.0
COLUMN_COEFF_MAX = 5.0

#==============================================================================

random.seed(42)

reflectors = np.zeros((GROUND_NUM_POINTS + COLUMN_NUM_POINTS, 4))

for i in range(GROUND_NUM_POINTS):
    reflectors[i, 0] = random.uniform(GROUND_X_MIN, GROUND_X_MAX)
    reflectors[i, 1] = random.uniform(GROUND_Y_MIN, GROUND_Y_MAX)
    reflectors[i, 2] = random.uniform(GROUND_Z_MIN, GROUND_Z_MAX)
    reflectors[i, 3] = random.uniform(GROUND_COEFF_MIN, GROUND_COEFF_MAX)

offset = GROUND_NUM_POINTS
for i in range(COLUMN_NUM_POINTS):
    r = abs(random.gauss(0.0, COLUMN_MAX_RADIUS / 2.0))
    theta = random.uniform(0.0, 2.0 * np.pi)
    z0 = (COLUMN_Z_MIN + COLUMN_Z_MAX) / 2.0
    z = abs(random.gauss(z0, 0.5 * (COLUMN_Z_MAX - COLUMN_Z_MIN) / 2.0))
    if z > COLUMN_Z_MAX: z = COLUMN_Z_MAX
    reflectors[offset + i, 0] = COLUMN_X_CENTER + r * np.cos(theta)
    reflectors[offset + i, 1] = COLUMN_Y_CENTER + r * np.sin(theta)
    reflectors[offset + i, 2] = z
    reflectors[offset + i, 3] = COLUMN_COEFF_MIN + (COLUMN_COEFF_MAX - COLUMN_COEFF_MIN) * random.random()

hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors-random_column_and_ground.h5", DATASET_NAME)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(reflectors[:, 0], reflectors[:, 1], reflectors[:, 2])
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
