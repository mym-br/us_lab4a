#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from util import arrayutil, hdf5util, decimation

DECIMATION_FILTER_TOLERANCE = 0.0001
DECIMATION_FILTER_TRANSITION_WIDTH = 0.3 # fraction of the original bandwidth
FS = 400e6
DECIMATION_OFFSET = 32
#DECIMATION_OFFSET = 33



ascan = hdf5util.read_to_ndarray(file_path='ref_pulse.h5', dataset_name='ascan')
ascan = np.hstack((ascan, np.zeros(500)))
ascan = np.hstack((ascan, np.ones(100)))
ascan = np.hstack((ascan, np.zeros(500)))

t = arrayutil.get_time_sequence(ascan, FS, offset=0)

hdf5util.write_ndarray(data=ascan, file_path='decimation_source.h5', dataset_name='v')

for decimation_factor in [2, 5, 10]:
    print('### decimation_factor: {}'.format(decimation_factor))
    decim_filter = decimation.downsampling_filter(decimation_factor=decimation_factor,
                                                  transition_width=DECIMATION_FILTER_TRANSITION_WIDTH,
                                                  tolerance=DECIMATION_FILTER_TOLERANCE,
                                                  plot=False)
    print('filter length: {}'.format(len(decim_filter)))

    ascan_d_offset, ascan_d, t_d = decimation.decimate(DECIMATION_OFFSET,
                                                       ascan,
                                                       decimation_factor,
                                                       decim_filter,
                                                       t)
    print('len(ascan) {}'.format(len(ascan)))
    print('len(ascan_d) {}'.format(len(ascan_d)))

    plt.figure()
    plt.plot(t, ascan, 'o-', label='orig')
    plt.plot(t_d, ascan_d, '.-', label='decimated ' + str(decimation_factor) + 'x')
    plt.legend(loc='upper right', labelspacing=0.2)
    plt.title('[t] decimated ' + str(decimation_factor) + 'x')
    plt.grid()

    plt.figure()
    plt.plot(np.arange(len(ascan)) + DECIMATION_OFFSET, ascan, 'o-', label='orig')
    plt.plot(decimation_factor * (np.arange(len(ascan_d)) + ascan_d_offset), ascan_d, '.-', label='decimated ' + str(decimation_factor) + 'x')
    plt.legend(loc='upper right', labelspacing=0.2)
    plt.title('[index] decimated ' + str(decimation_factor) + 'x')
    plt.grid()

    hdf5util.write_ndarray(data=ascan_d, file_path='decimated_' + str(decimation_factor) + 'x.h5', dataset_name='v')

plt.show()
