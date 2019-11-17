# This file is in the public domain.

import numpy as np
from util import hdf5util
from scipy.signal import resample

def waveform1(fc, sample_rate, num_periods):
    if num_periods is None:
        num_periods = 3.0
    t = np.arange(0.0, num_periods * sample_rate / fc) / sample_rate
    w = 2.0 * np.pi * fc
    return np.sin(w * t) * (0.5 - 0.5 * np.cos(w * (1.0 / num_periods) * t))

# Emeterio, J. L. S.
# Ullate, L. G.
# Diffraction impulse response of rectangular transducers.
# J. Acoust. Soc. Am., vol. 92, no. 2, pp. 651-662, 1992.
# DOI: 10.1121/1.403990
def _waveform2(fc, sample_rate, num_periods, K):
    t = np.arange(0.0, num_periods * sample_rate / fc) / sample_rate
    w = 2.0 * np.pi * fc
    x = t**3 * np.exp(-K * fc * t) * np.cos(w * t)
    return x / np.max(np.abs(x))
def waveform2a(fc, sample_rate, num_periods):
    if num_periods is None:
        num_periods = 3.25
    return _waveform2(fc, sample_rate, num_periods, 3.833)
def waveform2b(fc, sample_rate, num_periods):
    if num_periods is None:
        num_periods = 8.25
    return _waveform2(fc, sample_rate, num_periods, 1.437)

def hdf5_waveform(file_path, dataset_name, source_sample_rate, sample_rate):
    s = hdf5util.read_to_ndarray(file_path, dataset_name)
    s_r = resample(s, len(s) * round(sample_rate / source_sample_rate))
    return s_r
