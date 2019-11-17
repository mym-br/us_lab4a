#!/usr/bin/env python3
# This file is in the public domain.

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

t = arrayutil.get_time_sequence(len(ascan), FS)

hdf5util.write_ndarray(data=ascan, file_path='decimation_source.h5', dataset_name='v')

for decimation_factor in [2, 5, 10]:
    print('### decimation_factor: {}'.format(decimation_factor))
    decim_filter = decimation.downsampling_filter(decimation_factor=decimation_factor,
                                                  half_transition_width=DECIMATION_FILTER_TRANSITION_WIDTH,
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

    out_file_prefix = 'decimated_' + str(decimation_factor) + 'x'
    hdf5util.write_ndarray(data=ascan_d, file_path=out_file_prefix+'.h5', dataset_name='v')
    with open(out_file_prefix + '-offset.txt', 'w') as f:
        print(str(ascan_d_offset), file=f)

plt.show()
