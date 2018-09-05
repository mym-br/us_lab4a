#!/usr/bin/env python3
# This file is in the public domain.

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from util import hdf5util, imagegrid

#------------------------------------------------------------------------------
# Configurable parameters.

DATA_DIR = "../project/simulated_3d_sta-imasonic_0.5mhz_64elem/saved_acquisition/0000/"
SIGNAL_DATASET = "signal"

# (Hz)
CENTER_FREQ = 0.5e6
MAX_FREQ = 2.0 * CENTER_FREQ
# (m/s)
PROPAGATION_SPEED = 1500.0
# (Hz)
SAMPLING_FREQ = 4.0e6

BASE_ELEMENT = 16
NUM_ELEMENTS = 32
# First element: 0
MIN_TX_ELEM = 15
MAX_TX_ELEM = 15
# (m)
PITCH = 1.5e-3

USE_ENVELOPE = True
USE_DB_LEVELS = True
MIN_DB_LEVEL = -60.0
INVERT_X = False
INVERT_Z = True

# (m)
SECTORIAL_MIN_RADIUS = 0.2
SECTORIAL_MAX_RADIUS = 1.0
# One wavelength (at Nyquist rate) divided by this number.
SECTORIAL_RADIUS_STEP_DIV = 16.0
# (degree)
SECTORIAL_MIN_ANGLE = -45.0
SECTORIAL_MAX_ANGLE = 45.0
SECTORIAL_ANGLE_STEP = 0.25
# (m)
SECTORIAL_ORIGIN_X = 0.0
SECTORIAL_ORIGIN_Z = 0.0

# End of configurable parameters.
#------------------------------------------------------------------------------

nyquist_rate = 2.0 * MAX_FREQ
nyquist_lambda = PROPAGATION_SPEED / nyquist_rate;

x, z = imagegrid.gen_sectorial_grid_xz(SECTORIAL_MIN_RADIUS,
                                       SECTORIAL_MAX_RADIUS,
                                       nyquist_lambda / SECTORIAL_RADIUS_STEP_DIV,
                                       SECTORIAL_MIN_ANGLE,
                                       SECTORIAL_MAX_ANGLE,
                                       SECTORIAL_ANGLE_STEP,
                                       SECTORIAL_ORIGIN_X,
                                       SECTORIAL_ORIGIN_Z)
nz, nx = x.shape
if USE_ENVELOPE:
    image = np.zeros(x.shape, dtype=np.complex128)
else:
    image = np.zeros(x.shape)
print("Image shape: {}".format(image.shape))
delays = np.zeros((NUM_ELEMENTS, nz, nx))

# Calculate the one-way delays.
for elem in range(NUM_ELEMENTS):
    x_elem = (elem - 0.5 * (NUM_ELEMENTS - 1)) * PITCH
    dx = x - x_elem
    delays[elem, :, :] = (SAMPLING_FREQ / PROPAGATION_SPEED) * np.sqrt(dx * dx + z * z)

# Delay and sum.
for itx in range(MIN_TX_ELEM, MAX_TX_ELEM + 1):
    print("TX elem: {}".format(itx))
    file_path = DATA_DIR + "signal-base{:04d}-tx{:04d}.h5".format(BASE_ELEMENT, itx)
    print("Signal file: " + file_path)
    signals = hdf5util.read_to_ndarray(file_path, SIGNAL_DATASET)
    if USE_ENVELOPE:
        signals = hilbert(signals)
    tx_delays = delays[itx, :, :]
    signal_len = signals.shape[1]
    if NUM_ELEMENTS != signals.shape[0]:
        raise ValueError("Wrong number of receive elements: {} (should be {}).".format(signals.shape[0], NUM_ELEMENTS))
    print("Signal length: {}".format(signal_len))
    for irx in range(NUM_ELEMENTS):
        print("RX elem: {}".format(irx))
        rx_delays = delays[irx, :, :]
        total_delays = tx_delays + rx_delays
        index = np.around(total_delays).astype(int)
        index[index >= signal_len] = signal_len - 1
        image += signals[irx, :][index]

# Normalize.
image = np.abs(image)
coef = 1.0 / image.max()
image *= coef

if USE_DB_LEVELS:
    min_level = 10.0**(MIN_DB_LEVEL / 20.0)
    image[image < min_level] = min_level
    image = 20.0 * np.log10(image)
if INVERT_X:
    x = -x
if INVERT_Z:
    z = -z

plt.figure(figsize=(10, 7))
plt.pcolormesh(x, z, image)
plt.axis("equal")
plt.grid(True)
plt.xlabel("x (m)")
plt.ylabel("z (m)")
cbar = plt.colorbar()
if USE_DB_LEVELS:
    cbar.ax.set_ylabel('dB')
plt.autoscale()
plt.tight_layout(pad=0.5)
#plt.savefig("figure.png", dpi=300)

plt.show()
