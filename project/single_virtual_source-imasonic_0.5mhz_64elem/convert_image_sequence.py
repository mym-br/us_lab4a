#!/usr/bin/env python3

import sys
sys.path.append("../../python")
from util import hdf5util
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

#------------------------------------------------------------------------------
# Configurable parameters.

DATA_DIR = "./output/"

X_FILE = "/image_x.h5"
X_FILE_DATASET = "x"

Z_FILE = "/image_z.h5"
Z_FILE_DATASET = "z"

IMAGE_FILE         = "/image_value.h5"
IMAGE_FILE_DATASET = "value"
#IMAGE_FILE         = "/image_factor.h5"
#IMAGE_FILE_DATASET = "factor"
#IMAGE_FILE         = "/image_cf.h5"
#IMAGE_FILE_DATASET = "cf"

INVERT_X = False
INVERT_Z = False

TIME_FILE = "time.h5"
TIME_DATASET = "time"

USE_ENVELOPE = False

USE_DB_LEVELS = True
MIN_DB_LEVEL = -60.0

# End of configurable parameters.
#------------------------------------------------------------------------------

acq_number_format = "{:04}"

x = hdf5util.read_to_ndarray(DATA_DIR + acq_number_format.format(0) + X_FILE, X_FILE_DATASET)
if INVERT_X:
    x = -x

z = hdf5util.read_to_ndarray(DATA_DIR + acq_number_format.format(0) + Z_FILE, Z_FILE_DATASET)
if INVERT_Z:
    z = -z

def convert_image(acq_number):
    print("ACQ {}".format(acq_number))

    acq_number_str = acq_number_format.format(acq_number)

    image = hdf5util.read_to_ndarray(DATA_DIR + acq_number_str + IMAGE_FILE,
                                     IMAGE_FILE_DATASET)
    if USE_ENVELOPE:
        image = np.abs(hilbert(image, axis=1))
    else:
        image = np.abs(image)

    # Normalize.
    image *= 1.0 / image.max()

    if USE_DB_LEVELS:
        min_level = 10.0**(MIN_DB_LEVEL / 20.0)
        image[image < min_level] = min_level
        image = 20.0 * np.log10(image)

    plt.figure(figsize=(12, 9))
    plt.pcolormesh(x, z, image)
    plt.axis("equal")
    plt.xlabel("x (m)")
    plt.ylabel("z (m)")
    #cbar = plt.colorbar()
    #if USE_DB_LEVELS:
    #    cbar.ax.set_ylabel('dB')
    plt.autoscale()
    plt.tight_layout(pad=0.5)
    plt.savefig(DATA_DIR + IMAGE_FILE[0:-3] + "-" + acq_number_str + ".png", dpi=300)
    plt.close()

# Turn the interactive mode off.
plt.ioff()

time_list = hdf5util.read_to_ndarray(DATA_DIR + TIME_FILE, TIME_DATASET)

for i, t in enumerate(time_list.T):
    #if i > 10: break
    convert_image(i)
