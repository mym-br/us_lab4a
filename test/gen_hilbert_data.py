#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from util import arrayutil, hdf5util
from scipy.signal import hilbert

MIN_PADDING = 1024

x = hdf5util.read_to_ndarray(file_path='base48_tx15_rx15-calib_plane_40mm_avg1.h5', dataset_name='ascan')
x = x.flatten()
print('len(x):', len(x))

#size = np.maximum(arrayutil.get_next_power_of_two(len(x) + MIN_PADDING), 2 * MIN_PADDING)
size = np.maximum(arrayutil.get_next_fast_even_fft_size(len(x) + MIN_PADDING), 2 * MIN_PADDING)
print('size:', size)

xp = arrayutil.pad_with_zeros(x, size)

h = hilbert(xp)[:len(x)]
print('h.shape:', h.shape)
ah = np.abs(h)

plt.figure()
plt.plot(x)
plt.title('x')
plt.grid(True)

plt.figure()
plt.plot(h.real, 'r')
plt.plot(h.imag, 'b')
plt.title('h')
plt.grid(True)

plt.figure()
plt.plot(ah, 'r')
plt.plot(np.abs(x), 'b')
plt.title('ah')
plt.grid(True)

hdf5util.write_ndarray(data=x, file_path='hilbert_x.h5', dataset_name='v')
hdf5util.write_ndarray(data=ah, file_path='hilbert_ya.h5', dataset_name='v')
hdf5util.write_ndarray(data=h.real, file_path='hilbert_yr.h5', dataset_name='v')
hdf5util.write_ndarray(data=h.imag, file_path='hilbert_yi.h5', dataset_name='v')

plt.show()
