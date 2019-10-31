#!/usr/bin/env python3
# This file is in the public domain.

import numpy as np
import matplotlib.pyplot as plt
import util.hdf5util as h5

#------------------------------------------------------------------------------
# Configurable parameters.

DATA_DIR = ("output-radiation_pattern-transient-numeric/")

PATTERN_FILE = "pattern_theta_x.h5"
PATTERN_DATASET = "value"

ANGLE_FILE = "theta_x.h5"
ANGLE_DATASET = "value"

USE_DB_LEVELS = False
MIN_DB_LEVEL = -20.0

# End of configurable parameters.
#------------------------------------------------------------------------------

p = h5.read_to_ndarray(DATA_DIR + PATTERN_FILE, PATTERN_DATASET).flatten()
a = np.deg2rad(h5.read_to_ndarray(DATA_DIR + ANGLE_FILE, ANGLE_DATASET)).flatten()

if USE_DB_LEVELS:
    # Normalize.
    coef = 1.0 / p.max()
    p *= coef

    min_level = 10.0**(MIN_DB_LEVEL / 20.0)
    p[p < min_level] = min_level
    p = 20.0 * np.log10(p)

fig = plt.figure(figsize=(6, 6))
plt.polar(a, p)
plt.yticks(np.r_[np.array([0.0, -3.0]), np.arange(-6.0, MIN_DB_LEVEL, -6.0)[::-1]])
plt.ylim(MIN_DB_LEVEL, 0.0)

xt = np.arange(0.0, 360.0, 10.0)
xt_labels = xt.copy()
xt_labels[xt_labels > 180.0] -= 360.0
plt.xticks(np.deg2rad(xt), xt_labels.astype(int))

ax = fig.gca()
ax.set_rlabel_position(86.0)

plt.tight_layout(pad=0.5)
#plt.savefig("radiation_pattern.png", dpi=300)
plt.savefig("radiation_pattern.pdf")

plt.show()
