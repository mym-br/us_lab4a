#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import rfft
from pathlib import Path
import sys
from util.arrayutil import get_time_sequence
from util.signal import filter_spikes, sinc_lowpass
from util import hdf5util

PLOT_SIGNAL_INDEX = 16

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

def show_usage():
    print("Usage:")
    print("{} single <in_dir> <out_dir> <hdf5_file> <hdf5_dataset> <sampling_freq> <cutoff_freq>".format(sys.argv[0]))
    print("{} sequence <in_dir> <out_dir> <hdf5_dataset> <sampling_freq> <cutoff_freq>".format(sys.argv[0]))
    print("  sampling_freq: Hz")
    print("  cutoff_freq:   Hz")

def process_file(file_path):
    if in_dir == out_dir:
        raise ValueError("The input and output directories must be different.")
    signal_set = hdf5util.read_to_ndarray(in_dir + "/" + file_path, dataset_name)
    signal_set_out = np.zeros(signal_set.shape)
    for i in range(signal_set.shape[0]):
         s_out = filter_signal(signal_set[i, :], i == plot_signal_index)
         signal_set_out[i, :len(s_out)] = s_out

    out_path = out_dir + "/" + file_path
    op = Path(out_path)
    op.parent.mkdir(parents=True, exist_ok=True)

    hdf5util.write_ndarray(signal_set_out, out_path, dataset_name)



if __name__ == "__main__":

    if len(sys.argv) < 2:
        show_usage()
        sys.exit(1)

    option = sys.argv[1]
    if option == "single":
        if len(sys.argv) != 8:
            show_usage()
            sys.exit(1)
        in_dir       =       sys.argv[2]
        out_dir      =       sys.argv[3]
        file_path    =       sys.argv[4]
        dataset_name =       sys.argv[5]
        fs           = float(sys.argv[6])
        cutoff_freq  = float(sys.argv[7])
        plot_signal_index = PLOT_SIGNAL_INDEX

        dt = 1.0 / fs

        process_file(file_path)
    elif option == "sequence":
        if len(sys.argv) != 7:
            show_usage()
            sys.exit(1)
        in_dir       =       sys.argv[2]
        out_dir      =       sys.argv[3]
        dataset_name =       sys.argv[4]
        fs           = float(sys.argv[5])
        cutoff_freq  = float(sys.argv[6])
        plot_signal_index = -1

        dt = 1.0 / fs

        d = Path(in_dir)
        if not d.exists():
            raise FileNotFoundError("The directory {} does not exist.".format(in_dir))
        if not d.is_dir():
            raise NotADirectoryError("The path {} is not a directory.".format(in_dir))
        for child in sorted(d.iterdir()):
            if child.is_dir():
                for item in sorted(child.glob("*.h5")):
                    file_path = str(child.name) + "/" + str(item.name)
                    print(file_path)
                    process_file(file_path)
    else:
        show_usage()
        sys.exit(1)

    plt.show()
