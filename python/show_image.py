#!/usr/bin/env python3
# This file is in the public domain.

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import hilbert

#------------------------------------------------------------------------------
# Configurable parameters.

#DATA_DIR = "../../us_lab4a-processing/lab4a/sta-imasonic-0.5_mhz-64_elem/output/"
DATA_DIR = "../../us_lab4a-processing/lab4a/simulated_3d_sta_acquisition/output/"

X_FILE = "image_x.h5"
X_FILE_DATASET = "x"

Z_FILE = "image_z.h5"
Z_FILE_DATASET = "z"

IMAGE_FILE = "image_value.h5"
IMAGE_FILE_DATASET = "value"

USE_ENVELOPE = True

USE_DB_LEVELS = True
MIN_DB_LEVEL = -50.0

# End of configurable parameters.
#------------------------------------------------------------------------------

h5f = h5py.File(DATA_DIR + X_FILE, "r")
x = h5f[X_FILE_DATASET][:] # gets a ndarray from a Dataset using the [:]
h5f.close()

h5f = h5py.File(DATA_DIR + Z_FILE, "r")
z = h5f[Z_FILE_DATASET][:] # gets a ndarray from a Dataset using the [:]
h5f.close()

h5f = h5py.File(DATA_DIR + IMAGE_FILE, "r")
image = h5f[IMAGE_FILE_DATASET][:] # gets a ndarray from a Dataset using the [:]
h5f.close()

if USE_ENVELOPE:
    image = np.abs(hilbert(image, axis=1))
else:
    image = np.abs(image)

# Normalize.
coef = 1.0 / image.max()
image *= coef

if USE_DB_LEVELS:
    min_level = 10.0**(MIN_DB_LEVEL / 20.0)
    image[image < min_level] = min_level
    image = 20.0 * np.log10(image)

plt.pcolormesh(-x, z, image)
plt.axis("equal")
plt.grid(True)
plt.xlabel("x (m)")
plt.ylabel("z (m)")
plt.colorbar()
plt.autoscale()

plt.show()
