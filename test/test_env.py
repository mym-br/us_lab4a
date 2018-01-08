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

ascan_dataset = h5py.File('ascan.h5')['signal']
ascan = np.zeros(ascan_dataset.shape)
ascan_dataset.read_direct(ascan)

#plt.figure(1, figsize=(8, 3.5), dpi=85)
#plt.clf()
#plt.hold(True)
plt.plot(np.abs(ascan[15, :]))
plt.grid(True)
plt.hold(True)

ascan_env_dataset = h5py.File('ascan_env.h5')['signal']
ascan_env = np.zeros(ascan_env_dataset.shape)
ascan_env_dataset.read_direct(ascan_env)
plt.plot(ascan_env[15, :])

plt.show()
