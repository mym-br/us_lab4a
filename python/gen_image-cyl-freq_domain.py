#!/usr/bin/env python3
# This file is in the public domain.

# Only one transmit element.
# Fills with zero to interpolate.

import math
import numpy as np
import matplotlib.pyplot as plt
from util import hdf5util
from numpy.fft import fft2, ifft2

#------------------------------------------------------------------------------
# Configurable parameters.

DATA_DIR = "../project/simulated_3d_sta-imasonic_0.5mhz_64elem/saved_acquisition/0000/"
SIGNAL_DATASET = "signal"
REFLECTORS_FILE = "../project/simulated_3d_sta-imasonic_0.5mhz_64elem/reflectors0.h5"
REFLECTORS_DATASET = "reflectors"

# (m/s)
PROPAGATION_SPEED = 1500.0
# (Hz)
SAMPLING_FREQ = 4.0e6

BASE_ELEMENT = 16
NUM_ELEMENTS = 32
# First element: 0
TX_ELEM = 15
# (m)
PITCH = 1.5e-3

USE_DB_LEVELS = True
MIN_DB_LEVEL = -60.0
INVERT_X = False
INVERT_Z = True

W_OVERSAMP_FACTOR = 4
KV_OVERSAMP_FACTOR = 4

NX =   1024 # grid x size - must be power of two (?)
NZ = 2*1024 # grid z size - must be power of two (?)

# (m)
X_MIN = -0.5
X_MAX = 0.5
Z_MIN = 0.6
Z_MAX = 0.9
#X_MIN = -0.2
#X_MAX = 0.2
#Z_MIN = 0.05
#Z_MAX = 0.2

# End of configurable parameters.
#------------------------------------------------------------------------------

x_step = (X_MAX - X_MIN) / (NX - 1)
z_step = (Z_MAX - Z_MIN) / (NZ - 1)
two_pi = 2.0 * math.pi
halfWidth = 0.5 * (NUM_ELEMENTS - 1) * PITCH

reflectors = hdf5util.read_to_ndarray(REFLECTORS_FILE, REFLECTORS_DATASET)
if INVERT_X: reflectors[:, 0] = -reflectors[:, 0]
if INVERT_Z: reflectors[:, 2] = -reflectors[:, 2]

# Load the captured signals.
file_path = DATA_DIR + "signal-base{:04d}-tx{:04d}.h5".format(BASE_ELEMENT, TX_ELEM)
print("Signal file: " + file_path)
signals = hdf5util.read_to_ndarray(file_path, SIGNAL_DATASET)
signal_len = signals.shape[1]
if NUM_ELEMENTS != signals.shape[0]:
    raise ValueError("Wrong number of receive elements: {} (should be {}).".format(signals.shape[0], NUM_ELEMENTS))
print("Signal length: {}".format(signal_len))

min_sample = 0
max_sample = math.floor(SAMPLING_FREQ * 2.0 * Z_MAX / PROPAGATION_SPEED)
if max_sample >= signal_len:
    max_sample = signal_len - 1
num_valid_samples = max_sample - min_sample + 1
if num_valid_samples < 1:
    raise ValueError("Number of valid samples < 1.")

signals_os = np.zeros((NUM_ELEMENTS * KV_OVERSAMP_FACTOR, signal_len * W_OVERSAMP_FACTOR))
for elem in range(NUM_ELEMENTS):
    signals_os[elem, :signal_len] = signals[elem, :]

plt.figure()
plt.pcolormesh(signals_os)
plt.title("signals_os")

# kv (receiver), angular velocity w
s_kvw = fft2(signals_os)

#plt.figure()
#plt.pcolormesh(s_kvw.real)
#plt.title("s_kvw.real")
#
#plt.figure()
#plt.pcolormesh(s_kvw.imag)
#plt.title("s_kvw.imag")

plt.figure()
plt.pcolormesh(abs(s_kvw))
plt.title("abs(s_kvw)")

#------------------------------------------------------------------------------

s_kxkz = np.zeros((NX, NZ), dtype=complex)
s_kxkz_mask = np.zeros(s_kxkz.shape)
s_kxkz_transl = np.zeros(s_kxkz.shape, dtype=complex)

w_period = two_pi * SAMPLING_FREQ
kv_period = two_pi / PITCH

half_nx = NX // 2
half_nz = NZ // 2

inv_x_len = 1.0 / (x_step * NX)
inv_z_len = 1.0 / (z_step * NZ)

two_pi_inv_x_len = two_pi * inv_x_len
two_pi_inv_z_len = two_pi * inv_z_len

ku = 0.0
ku2 = ku * ku

ikx_top_left = 0
ikz_top_left = 0
ikz_top_right = 0

#max_abs_k = 0.0
#max_abs_dk = 0.0
min_ikv = NUM_ELEMENTS * KV_OVERSAMP_FACTOR + 1
max_ikv = -1

for ikx in range(NX):
    if ikx > half_nx:
        kx = two_pi_inv_x_len * (ikx - NX)
    else:
        kx = two_pi_inv_x_len * ikx
    kx2 = kx * kx

    kv = kx - ku
    kv2 = kv * kv

    if abs(kv) >= kv_period / 2.0:
        print("CONTINUE 2 ikx={}".format(ikx))
        # In s_kxkz_mask, causes the horizontal central empty band.
        continue
    if kv < 0.0:
        ikv = (NUM_ELEMENTS * KV_OVERSAMP_FACTOR) * (kv + kv_period) / kv_period # float
    else:
        ikv = (NUM_ELEMENTS * KV_OVERSAMP_FACTOR) * kv / kv_period # float

    if ikv < min_ikv: min_ikv = ikv
    if ikv > max_ikv: max_ikv = ikv

    ikv_base = math.floor(ikv)
    ikv_factor = ikv - ikv_base

    for ikz in range(NZ // 2):
        #if ikz > half_nz:
        #    kz = two_pi_inv_z_len * (ikz - NZ)
        #else:
        kz = two_pi_inv_z_len * ikz

        kz2 = kz * kz

        if kz == 0.0:
            print("CONTINUE 3 ikx={}".format(ikx))
            continue
        else:
            if ku == 0.0:
                k = (kx2 + kz2) / (2.0 * kz)
            else:
                k = math.sqrt((kz2 + 2.0 * (ku2 + kv2)) * kz2 + ku2 * ku2 + kv2 * kv2 - 2.0 * ku2 * kv2) / (2.0 * kz)

#            abs_dk = abs(k - k2)
#            if abs_dk > max_abs_dk:
#                max_abs_dk = abs_dk
#            abs_k = abs(k)
#            if abs_k > max_abs_k:
#                max_abs_k = abs_k

        w = k * PROPAGATION_SPEED
        if abs(w) >= w_period / 2.0:
            print("CONTINUE 5 ikx={}".format(ikx))
            continue
        if w < 0.0:
            raise ValueError("Negative frequency.")
            #iw = (signal_len * W_OVERSAMP_FACTOR) * (w + w_period) / w_period # float
        else:
            iw = (signal_len * W_OVERSAMP_FACTOR) * w / w_period # float

        iw_base = math.floor(iw)
        iw_factor = iw - iw_base

        if iw_base >= 0 and iw_base < signal_len * W_OVERSAMP_FACTOR - 1 \
                and ikv_base >= 0 and ikv_base < NUM_ELEMENTS * KV_OVERSAMP_FACTOR - 1:

            if ikx < NX // 2:
                if ikx > ikx_top_left:
                    ikx_top_left = ikx
                    ikz_top_left = ikz
                    ikz_top_right = ikz
                elif ikx == ikx_top_left:
                    if ikz < ikz_top_left:
                        ikz_top_left = ikz
                    if ikz > ikz_top_right:
                        ikz_top_right = ikz

            k2 = k * k
#            if k2.imag != 0.0:
#                raise ValueError("k2.imag != 0.0")

            # Bilinear interpolation.
            s_kvw_interp = \
                (1.0 - ikv_factor) * ((1.0 - iw_factor) * s_kvw[ikv_base    , iw_base] + iw_factor * s_kvw[ikv_base    , iw_base + 1]) + \
                ikv_factor         * ((1.0 - iw_factor) * s_kvw[ikv_base + 1, iw_base] + iw_factor * s_kvw[ikv_base + 1, iw_base + 1])

            if k2 >= ku2 and k2 >= kv2:
                s_kxkz[ikx, ikz] = -math.sqrt(k2 - ku2) * math.sqrt(k2 - kv2) * s_kvw_interp

                s_kxkz_mask[ikx, ikz] = 1.0
                # Translation compensation.
                s_kxkz_transl[ikx, ikz] = np.exp(-1j * (-X_MIN * kx - Z_MIN * kz))
            else:
                raise ValueError("k2 < ku2 or k2 < kv2")

#print("max_abs_dk={}".format(max_abs_dk))
#print("max_abs_k={}".format(max_abs_k))
print("min_ikv={}".format(min_ikv))
print("max_ikv={}".format(max_ikv))

# The original zero is the center of the element zero.
# Subtract halfWidth to move the zero to the center of the array.
x_list = np.linspace(X_MIN, X_MAX, NX) - halfWidth
z_list = np.linspace(Z_MIN, Z_MAX, NZ)
if INVERT_X: x_list = -x_list
if INVERT_Z: z_list = -z_list
gx, gz = np.meshgrid(x_list, z_list, indexing='ij')

#plt.figure()
#plt.pcolormesh(s_kxkz.real)
#plt.title("s_kxkz.real")
#
#plt.figure()
#plt.pcolormesh(s_kxkz.imag)
#plt.title("s_kxkz.imag")

plt.figure()
plt.pcolormesh(np.abs(s_kxkz))
plt.title("abs(s_kxkz)")

plt.figure()
plt.pcolormesh(s_kxkz_mask)
plt.title("s_kxkz_mask")

s_kxkz_mask2 = s_kxkz_mask.copy()
s_kxkz_mask2[:, 0:ikz_top_left-1] = 0.0
s_kxkz_mask2[:, ikz_top_right+1:] = 0.0
plt.figure()
plt.pcolormesh(s_kxkz_mask2)
plt.title("s_kxkz_mask 2")

s_kxkz2 = s_kxkz.copy()
s_kxkz2[:, 0:ikz_top_left-1] = 0.0
s_kxkz2[:, ikz_top_right+1:] = 0.0
# Apply translation.
s_kxkz2 *= s_kxkz_transl
image = ifft2(s_kxkz2)

#plt.figure()
#plt.pcolormesh(gx, gz, image.real)
#plt.title("image.real")
#
#plt.figure()
#plt.pcolormesh(gx, gz, image.imag)
#plt.title("image.imag")

image_abs = abs(image)
# Normalize.
coef = 1.0 / image_abs.max()
image_abs *= coef
if USE_DB_LEVELS:
    min_level = 10.0**(MIN_DB_LEVEL / 20.0)
    image_abs[image_abs < min_level] = min_level
    image_abs = 20.0 * np.log10(image_abs)

plt.figure()
plt.pcolormesh(gx, gz, image_abs)
plt.plot(reflectors[:, 0], reflectors[:, 2], 'wx')
plt.title("image")

print("w_period: {}".format(w_period))
print("kv_period: {}".format(kv_period))
print("x_step: {}".format(x_step))
print("z_step: {}".format(z_step))

plt.show()
