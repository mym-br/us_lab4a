#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import rfft
import sys
from util.arrayutil import get_time_sequence
from util.signal import filter_spikes, sinc_lowpass
from util import hdf5util

def plot_spectrum(signal, title):
    f = rfft(signal)
    n = len(signal)
    i = np.arange(n)
    freq = ((i / n) * fs)[:len(f)]

    plt.figure()
    plt.title(title + ' - abs')
    plt.plot(freq, abs(f))
    plt.xlabel('freq (Hz)')
    plt.grid(True)

    plt.figure()
    plt.title(title + ' - real')
    plt.plot(freq, f.real)
    plt.xlabel('freq (Hz)')
    plt.grid(True)

    plt.figure()
    plt.title(title + ' - imag')
    plt.plot(freq, f.imag)
    plt.xlabel('freq (Hz)')
    plt.grid(True)

def filter_signal(signal, plot_signals):
    s2 = filter_spikes(signal)
    s_filtered = sinc_lowpass(s2, fs, cutoff_freq)
    if plot_signals:
        plt.figure()
        plt.title("Comparison of signals")
        plt.plot(get_time_sequence(signal    , fs, 0.0), signal    , "r-o", label="raw")
        plt.plot(get_time_sequence(s2        , fs, 0.0), s2        , "g--", label="spikes filtered")
        plt.plot(get_time_sequence(s_filtered, fs, 0.0), s_filtered, "b--", label="filtered")
        plt.legend(loc='upper right')
        plt.xlabel('t (s)')
        plt.grid(True)

        plot_spectrum(signal    , "Spetrum - raw signal")
        plot_spectrum(s_filtered, "Spectrum - filtered signal")
    return s_filtered



if __name__ == "__main__":

    if len(sys.argv) != 5:
        print("Usage: {} hdf5_file hdf5_dataset sampling_freq cutoff_freq".format(sys.argv[0]))
        print("  sampling_freq: Hz")
        print("  cutoff_freq:   Hz")
        sys.exit(1)

    file_path    =       sys.argv[1]
    dataset_name =       sys.argv[2]
    fs           = float(sys.argv[3])
    cutoff_freq  = float(sys.argv[4])

    dt = 1.0 / fs

    signal_set = hdf5util.read_to_ndarray(file_path, dataset_name)
    signal_set_out = np.zeros(signal_set.shape)
    for i in range(signal_set.shape[0]):
         s_out = filter_signal(signal_set[i, :], i == 16)
         signal_set_out[i, :len(s_out)] = s_out

    hdf5util.write_ndarray(signal_set_out, "filtered_" + file_path, dataset_name)

    plt.show()
