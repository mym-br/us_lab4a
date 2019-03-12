#!/usr/bin/env python
# This file is in the public domain.

import h5py
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.signal

#rc("text", usetex=True)
#rc("font", family="serif")
#rc("font", size="15")
rc("font", size="14")
#matplotlib.rcParams['lines.linewidth'] = 1

plt.close('all')

plt.figure(1)
b_dataset = h5py.File('filter_b.h5')['b']
b = np.zeros(b_dataset.shape)
b_dataset.read_direct(b)
plt.plot(b[0, :])
plt.title('b')
print 'b.shape = ', b.shape
#plt.figure(1, figsize=(8, 3.5), dpi=85)
#plt.clf()
#plt.hold(True)

plt.figure(2)
x_dataset = h5py.File('filter_x.h5')['x']
x = np.zeros(x_dataset.shape)
x_dataset.read_direct(x)
plt.plot(x[0, :])
plt.title('x')
#plt.plot(b)
#plt.grid(True)
#plt.hold(True)

plt.figure(3)
y = sp.signal.convolve(x, b)
print 'y.shape = ', y.shape
plt.plot(y[0, :])
plt.title('y')
f = h5py.File('filter_y.h5', 'w')
dset = f.create_dataset('y', data=y, chunks=y.shape, compression='gzip', compression_opts=9)
f.close()

plt.figure(4)
y2 = sp.signal.convolve(x, y)
print 'y2.shape = ', y2.shape
plt.plot(y2[0, :])
plt.title('y2')
f = h5py.File('filter_y2.h5', 'w')
dset = f.create_dataset('y', data=y2, chunks=y2.shape, compression='gzip', compression_opts=9)
f.close()

plt.show()
