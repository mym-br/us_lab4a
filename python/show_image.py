#!/usr/bin/env python3
# This file is in the public domain.

import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import hilbert
#import matplotlib as mpl
#mpl.rc('image', cmap='inferno')

#------------------------------------------------------------------------------
# Configurable parameters.

DATA_DIR = "../project/simulated_3d_sta-imasonic_0.5mhz_64elem/output/"
#DATA_DIR = "../project/sta_network_sync-imasonic_0.5mhz_64elem/output/0114/"

X_FILE = "image_x.h5"
X_FILE_DATASET = "x"

Z_FILE = "image_z.h5"
Z_FILE_DATASET = "z"

IMAGE_FILE = "image_value.h5"
IMAGE_FILE_DATASET = "value"
#IMAGE_FILE = "image_cf.h5"
#IMAGE_FILE_DATASET = "cf"

USE_ENVELOPE = False

USE_DB_LEVELS = True
MIN_DB_LEVEL = -60.0

INVERT_X = False
INVERT_Z = True

# End of configurable parameters.
#------------------------------------------------------------------------------

h5f = h5py.File(DATA_DIR + X_FILE, "r")
x = h5f[X_FILE_DATASET][:] # gets a ndarray from a Dataset using the [:]
h5f.close()
if INVERT_X:
    x = -x

h5f = h5py.File(DATA_DIR + Z_FILE, "r")
z = h5f[Z_FILE_DATASET][:] # gets a ndarray from a Dataset using the [:]
h5f.close()
if INVERT_Z:
    z = -z

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

plt.figure(figsize=(10, 7))
plt.pcolormesh(-x, z, image)
plt.axis("equal")
#plt.grid(True)
plt.xlabel("x (m)")
plt.ylabel("z (m)")
cbar = plt.colorbar()
if USE_DB_LEVELS:
    cbar.ax.set_ylabel('dB')
plt.autoscale()
plt.tight_layout(pad=0.5)
plt.savefig("figure.png", dpi=300)

plt.show()
