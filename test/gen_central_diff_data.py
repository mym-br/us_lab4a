#!/usr/bin/env python
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from util import arrayutil, hdf5util

FS = 40e6

x = hdf5util.read_to_ndarray(file_path='base48_tx15_rx15-calib_plane_40mm_avg1.h5', dataset_name='ascan')
x = x.flatten()
#x = np.array([0., 0, 1, 1, 1, 1])
print x.shape

t = arrayutil.get_time_sequence(x, FS, 0)
t_d, x_d = arrayutil.central_diff_with_time(t, x, FS)

print t_d.shape
print t[0], t_d[0]
print x_d.shape

dx = arrayutil.central_diff(x, FS)
print dx.shape

print 'max error:', np.abs(x_d - dx).max()

plt.figure()
plt.plot(x)

plt.figure()
plt.plot(dx)

hdf5util.write_ndarray(data=x, file_path='central_diff_x.h5', dataset_name='v')
hdf5util.write_ndarray(data=dx, file_path='central_diff_y.h5', dataset_name='v')

plt.show()
