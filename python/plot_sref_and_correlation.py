#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from util import arrayutil, hdf5util
from numpy import correlate

from util.plotutil import PlotConfig
pcfg = PlotConfig('fig/plot_sref_and_correlation')
pcfg.set_tick_labels_padding(6.0)

ASCAN_FILE_NAME = '../data/ref_pulse-40MHz.h5'
ASCAN_DATASET_NAME = 'ascan'
FS = 40.0e6

PAD_TIME = 2.0e-6

#------------------------------------------------------------------------------

#fs = FS * UPSAMPLING_FACTOR
fs = FS

pad_samples = int(PAD_TIME * fs)

ascan = hdf5util.read_to_ndarray(ASCAN_FILE_NAME, ASCAN_DATASET_NAME)
t = arrayutil.get_time_sequence(len(ascan), FS, PAD_TIME)

ascan2 = np.r_[np.zeros(pad_samples), ascan, np.zeros(pad_samples)]
t2 = arrayutil.get_time_sequence(len(ascan2), fs, 0)

xc = correlate(ascan, ascan2, 'valid')
t_xc = arrayutil.get_time_sequence(len(xc), fs, 0)

plt.figure(figsize=(6, 3))
plt.plot(t2[:-pad_samples//2], ascan2[:-pad_samples//2] / np.abs(ascan2[:-pad_samples//2]).max(), 'k')
plt.plot(t, ascan / np.abs(ascan).max() - 2.0, 'k')
plt.plot(t_xc, xc / np.abs(xc).max() - 4.0, 'k')
plt.plot([PAD_TIME, PAD_TIME], [-5.0, 1.0], 'k--')
#plt.grid(True)
if pcfg.t1():
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
elif pcfg.t3():
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/ascan.svg')

plt.show()
