#!/usr/bin/env python3
# This file is in the public domain.

import sys
from util import hdf5util
import matplotlib.pyplot as plt
import numpy as np

DATASET_NAME = "apod"

#==============================================================================

def show_usage():
    print("Usage:")
    print("{} <window_size>".format(sys.argv[0]))

if __name__ == "__main__":

    if len(sys.argv) != 2:
        show_usage()
        sys.exit(1)

    window_size = int(sys.argv[1])
    if window_size <= 0:
        window_size = 1
    elif window_size > 4096:
        window_size = 4096

    apod = np.ones(window_size)

    hdf5util.write_ndarray(apod, "apod_1d_rectangular-" + str(window_size) + ".h5", DATASET_NAME)

    plt.stem(apod)
    plt.show()
