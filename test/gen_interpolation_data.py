#!/usr/bin/env python
# -*- coding: iso-8859-1 -*-

import matplotlib.pyplot as plt
import numpy as np
from util import arrayutil, hdf5util, interpolation

UPSAMP_FILTER_TOLERANCE = 0.0001
UPSAMP_FILTER_TRANSITION_WIDTH = 0.45 # fraction of the original bandwidth
FS = 40e6



ascan = hdf5util.read_to_ndarray(file_path='base48_tx15_rx15-calib_plane_40mm_avg1.h5', dataset_name='ascan')
ascan = ascan.T[2000:3000].flatten()
t = arrayutil.get_time_sequence(ascan, FS, offset=0)

hdf5util.write_ndarray(data=ascan, file_path='interp_source.h5', dataset_name='v')

for upsamp_factor in [2, 8, 64]:
    print '### upsamp_factor:', upsamp_factor
    upsamp_filter = interpolation.upsampling_filter(upsamp_factor=upsamp_factor,
                                                    transition_width=UPSAMP_FILTER_TRANSITION_WIDTH,
                                                    tolerance=UPSAMP_FILTER_TOLERANCE,
                                                    plot=False)
    filter_len = len(upsamp_filter)
    print 'filter length:', filter_len
    ascan_r, t_r = interpolation.interpolate(ascan,
                                             upsamp_factor,
                                             upsamp_filter,
                                             t,
                                             FS)
    offset = (len(ascan_r) - (len(ascan) * upsamp_factor)) / 2
    print 'offset:', offset
    ascan_r = ascan_r[offset:-offset]
    t_r = t_r[offset:-offset]
    print 'len(ascan)', len(ascan)
    print 'len(ascan_r)', len(ascan_r)

    max_error = np.abs(ascan - ascan_r[::upsamp_factor]).max()
    print 'max_error:', max_error

    plt.figure()
    plt.plot(t, ascan, label='orig')
    plt.plot(t_r, ascan_r, label='interpolated ' + str(upsamp_factor) + 'x')
    plt.legend(loc='upper right', labelspacing=0.2)
    plt.grid()

    plt.figure()
    plt.plot(np.arange(len(ascan_r)), ascan_r)
    plt.title('interpolated ' + str(upsamp_factor) + 'x')
    plt.grid()

    hdf5util.write_ndarray(data=ascan_r, file_path='interp_' + str(upsamp_factor) + 'x.h5', dataset_name='v')

plt.show()
