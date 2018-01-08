#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import h5py
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

#rc("text", usetex=True)
#rc("font", family="serif")
#rc("font", size="15")
rc("font", size="14")
#matplotlib.rcParams['lines.linewidth'] = 1

plt.close('all')

interp_x_dataset = h5py.File('interp_x.h5')['x']
interp_x = np.zeros(interp_x_dataset.shape)
interp_x_dataset.read_direct(interp_x)

#plt.figure(1, figsize=(8, 3.5), dpi=85)
#plt.clf()
s = interp_x.shape

plt.plot(np.arange(0, s[1] * 4, 4), interp_x[0, :])
plt.grid(True)
plt.hold(True)

interp_y_dataset = h5py.File('interp_y.h5')['y']
interp_y = np.zeros(interp_y_dataset.shape)
interp_y_dataset.read_direct(interp_y)
plt.plot(np.arange(0, s[1] * 4), interp_y[0, :])

plt.show()
