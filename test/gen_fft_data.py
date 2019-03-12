#!/usr/bin/env python
# This file is in the public domain.

import matplotlib.pyplot as plt
#import numpy as np
from util import hdf5util
from scipy.fftpack import fft

x = hdf5util.read_to_ndarray(file_path='base48_tx15_rx15-calib_plane_40mm_avg1.h5', dataset_name='ascan')
x1 = x.flatten()
x1_f = fft(x1)

x2 = x1[:2048]
x2_f = fft(x2)

x3 = x1[:3571] # prime
x3_f = fft(x3)



plt.figure()
plt.plot(x1)
plt.title('x')

plt.figure()
plt.plot(x1_f.real, 'r')
plt.plot(x1_f.imag, 'b')
plt.title('x1_f')

plt.figure()
plt.plot(x2_f.real, 'r')
plt.plot(x2_f.imag, 'b')
plt.title('x2_f')

plt.figure()
plt.plot(x3_f.real, 'r')
plt.plot(x3_f.imag, 'b')
plt.title('x3_f')

hdf5util.write_ndarray(data=x1, file_path='fft_x_4000.h5', dataset_name='v')
hdf5util.write_ndarray(data=x2, file_path='fft_x_2048.h5', dataset_name='v')
hdf5util.write_ndarray(data=x3, file_path='fft_x_3571.h5', dataset_name='v')

hdf5util.write_ndarray(data=x1_f.real, file_path='fft_yr_4000.h5', dataset_name='v')
hdf5util.write_ndarray(data=x1_f.imag, file_path='fft_yi_4000.h5', dataset_name='v')

hdf5util.write_ndarray(data=x2_f.real, file_path='fft_yr_2048.h5', dataset_name='v')
hdf5util.write_ndarray(data=x2_f.imag, file_path='fft_yi_2048.h5', dataset_name='v')

hdf5util.write_ndarray(data=x3_f.real, file_path='fft_yr_3571.h5', dataset_name='v')
hdf5util.write_ndarray(data=x3_f.imag, file_path='fft_yi_3571.h5', dataset_name='v')

plt.show()
