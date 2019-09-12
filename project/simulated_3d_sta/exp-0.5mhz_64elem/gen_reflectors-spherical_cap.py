#!/usr/bin/env python3
# This file is in the public domain.

import sys
sys.path.append("../../../python")
import numpy as np
from util import hdf5util
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

PROJECT_DIR = "./"
DATASET_NAME = "reflectors"

# m
X_CENTER = 0.0
Y_CENTER = 0.0
Z_CENTER = 0.7
RADIUS = 5.0e-3
# degree (< 90)
ANGLE = 20.0
# m
GRID_STEP = 1.5e-3/4.0

Z_COEF = -1.0

#==============================================================================

wx_max = np.deg2rad(ANGLE)
wx_list = np.arange(0.0, wx_max, GRID_STEP / RADIUS)
wx_list = np.r_[-wx_list[:0:-1], wx_list]

point_list = []
z_min = RADIUS * np.cos(wx_max)
for wx in wx_list:
    z_max = RADIUS * np.cos(wx)
    wy_max = np.arccos(z_min / z_max)
    wy_list = np.arange(0.0, wy_max, GRID_STEP / z_max)
    wy_list = np.r_[-wy_list[:0:-1], wy_list]
    y = Y_CENTER + RADIUS * np.sin(wx)
    for wy in wy_list:
        point_list.append([X_CENTER + z_max * np.sin(wy),
                           y,
                           Z_CENTER + (z_max * np.cos(wy)) * Z_COEF,
                           1.0])

reflectors = np.array(point_list)

hdf5util.write_ndarray(reflectors, PROJECT_DIR + "reflectors-spherical_cap.h5", DATASET_NAME)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(reflectors[:, 0], reflectors[:, 1], reflectors[:, 2])
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
ax.set_zlabel('z (m)')
