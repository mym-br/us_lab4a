#!/usr/bin/env python3
# This file is in the public domain.

import sys
sys.path.append("../../../../python")
from util import hdf5util
import matplotlib.pyplot as plt

from glob import glob

#------------------------------------------------------------------------------
# Configurable parameters.

DATA_DIR = "./"

X_FILE = "image_x.h5"
X_FILE_DATASET = "x"

Z_FILE = "image_z.h5"
Z_FILE_DATASET = "z"

IMAGE_FILE_PREFIX = "image_value-"
IMAGE_FILE_DATASET = "value"

INVERT_X = False
INVERT_Z = False

# End of configurable parameters.
#------------------------------------------------------------------------------

x = hdf5util.read_to_ndarray(DATA_DIR + X_FILE, X_FILE_DATASET)
if INVERT_X:
    x = -x

z = hdf5util.read_to_ndarray(DATA_DIR + Z_FILE, Z_FILE_DATASET)
if INVERT_Z:
    z = -z

for f in glob(DATA_DIR + IMAGE_FILE_PREFIX + "*.h5"):
    print(f)

    image = hdf5util.read_to_ndarray(f, IMAGE_FILE_DATASET)

    plt.figure(figsize=(12, 9))
    plt.pcolormesh(x, z, image, cmap="gray")
    plt.axis("equal")
    plt.xlabel("x (m)")
    plt.ylabel("z (m)")
    plt.autoscale()
    plt.tight_layout(pad=0.5)
    plt.clim(-1.0, 1.0)
    plt.savefig(f[:-2] + "png", dpi=300)
    plt.close()
