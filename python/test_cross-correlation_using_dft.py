#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft, ifft
from util import arrayutil, hdf5util



PULSE_FILE_PATH = '../data/ref_pulse.h5'
PULSE_DATASET_NAME = 'ascan'

PADDING_SIZE_BEFORE = 4444
PADDING_SIZE_AFTER = 3333

#-------------------------------------------------------------------------------

pulse = hdf5util.read_to_ndarray(file_path=PULSE_FILE_PATH,
                                 dataset_name=PULSE_DATASET_NAME)
print('len(pulse):', len(pulse))
plt.figure()
plt.plot(pulse)
plt.grid()

signal = np.r_[np.zeros((PADDING_SIZE_BEFORE,)), pulse, np.zeros((PADDING_SIZE_AFTER,))]
plt.figure()
plt.plot(signal)
plt.grid()

rev_pulse = pulse[::-1] # reverse
plt.figure()
plt.plot(rev_pulse)
plt.grid()

n = arrayutil.get_next_power_of_two(max(len(signal), len(rev_pulse)))

if n > len(signal):
    signal = arrayutil.pad_with_zeros(signal, n)
if n > len(rev_pulse):
    rev_pulse = arrayutil.pad_with_zeros(rev_pulse, n)
plt.figure()
plt.plot(signal)
plt.plot(rev_pulse)
plt.grid()

s_f = fft(signal)
p_f = fft(rev_pulse)

cc = ifft(s_f * p_f)
print('max(abs(imag(cc))) / max(abs(cc))', np.abs(cc.imag).max() / np.abs(cc).max())
cc = cc.real
cc_peak_index = np.argmax(np.abs(cc))
print('cc peak index:', cc_peak_index)
plt.figure()
plt.plot(cc)
plt.grid()
plt.title('Cross-correlation')

print('echo position:', cc_peak_index - (len(pulse) - 1))

plt.show()
