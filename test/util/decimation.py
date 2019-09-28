# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import fftconvolve, freqz
from scipy.signal import kaiser, kaiserord



# transition_width:
#     half the total transition width
#     1.0 -> pi radian / sample at the destination sampling rate
def downsampling_filter(decimation_factor, transition_width, tolerance=0.1, plot=False):
    if (decimation_factor == 1): return [1.0]
    if transition_width > 1.0:
        raise ValueError('transition_width > 1.0')

    attenuation = -20.0 * np.log10(tolerance)

    w_size, beta = kaiserord(attenuation,
                             width=(transition_width * 2.0 / decimation_factor))
    print('window size: {} beta: {}'.format(w_size, beta))

    w_size = int(np.ceil((w_size - 1) / (2.0 * decimation_factor)) *
                 2.0 * decimation_factor + 1.0)
    print('new window size: {}'.format(w_size))

    w = kaiser(w_size, beta)
    num_periods = (w_size - 1) / (2.0 * decimation_factor)
    x = np.linspace(-num_periods, num_periods, w_size);
    s = np.sinc(x) / decimation_factor
    coef = s * w

    if plot:
        plt.figure(figsize=(12, 6), dpi=85)
        plt.stem(np.arange(0, len(w)), w)
        plt.title('Kaiser window size=' + str(w_size) + ' beta=' + str(beta))
        plt.grid(True)

        plt.figure(figsize=(12, 6), dpi=85)
        plt.stem(x, s)
        plt.title('Sinc period=' + str(decimation_factor))
        plt.grid(True)

        plt.figure(figsize=(12, 6), dpi=85)
        plt.stem(x, coef)
        plt.title('Windowed sinc')
        plt.grid(True)

        freq, h = freqz(coef, worN=2048)
        ft = 0.5 / decimation_factor
        plt.figure(figsize=(12, 6), dpi=85)
        cf = 1.0 / decimation_factor
        plt.plot(freq / (2.0 * np.pi), abs(h) / decimation_factor, 'k',
                 [0.0, ft], [cf + tolerance, cf + tolerance], 'r',
                 [0.0, ft], [cf - tolerance, cf - tolerance], 'r',
                 [ft, 0.5], [tolerance, tolerance], 'r',
                 [ft, ft], [0.0, cf], 'r')
        plt.title('Filter freq. response')
        plt.grid(True)

    return coef



def decimate(signal_offset, signal, decimation_factor, fir_lp_filter, t):
    signal_filtered = fftconvolve(signal, fir_lp_filter, mode="same")

    m = signal_offset % decimation_factor
    if m == 0:
        first_index = 0
    else:
        first_index = decimation_factor - m
    signal_out = signal_filtered[first_index::decimation_factor]
    t_out = t[first_index::decimation_factor]
    signal_offset_out = int((signal_offset + first_index) // decimation_factor)

    return signal_offset_out, signal_out, t_out
