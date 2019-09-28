# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import fftconvolve, freqz
from scipy.signal import kaiser, kaiserord



# transition_width:
#     half the total transition width
#     1.0 -> pi radian / sample at the original sampling rate
def upsampling_filter(upsamp_factor, transition_width, tolerance=0.1, plot=False):
    if (upsamp_factor == 1): return [1.0]
    if transition_width > 1.0:
        raise ValueError('transition_width > 1.0')

    attenuation = -20.0 * np.log10(tolerance)
    #print 'attenuation', attenuation

    #beta = kaiser_beta(FILTER_ATTENUATION)

    w_size, beta = kaiserord(attenuation,
                             width=(transition_width * 2.0 / upsamp_factor))
    #print 'window size:', w_size, 'beta:', beta

    w_size = int(np.ceil((w_size - 1) / (2.0 * upsamp_factor)) *
                 2.0 * upsamp_factor + 1.0)
    #print 'new window size:', w_size

    w = kaiser(w_size, beta)
    num_periods = (w_size - 1) / (2.0 * upsamp_factor)
    x = np.linspace(-num_periods, num_periods, w_size);
    s = np.sinc(x)
    coef = s * w

    if plot:
        plt.figure(figsize=(12, 6), dpi=85)
        plt.stem(np.arange(0, len(w)), w)
        plt.title('Kaiser window size=' + str(w_size) + ' beta=' + str(beta))

        plt.figure(figsize=(12, 6), dpi=85)
        plt.stem(x, s)
        plt.title('Sinc period=' + str(upsamp_factor))
        plt.grid(True)

        plt.figure(figsize=(12, 6), dpi=85)
        plt.stem(x, coef)
        plt.title('Windowed sinc')
        plt.grid(True)

        freq, h = freqz(coef, worN=2048)
        ft = 0.5 / upsamp_factor
        plt.figure(figsize=(12, 6), dpi=85)
        plt.plot(freq / (2.0 * np.pi), abs(h) / upsamp_factor, 'k',
                 [0.0, ft], [1.0 + tolerance, 1.0 + tolerance], 'r',
                 [0.0, ft], [1.0 - tolerance, 1.0 - tolerance], 'r',
                 [ft, 0.5], [tolerance, tolerance], 'r',
                 [ft, ft], [0.0, 1.0], 'r')
        plt.title('Filter freq. response')
        plt.grid(True)

    return coef



def interpolate(signal, upsamp_factor, fir_lp_filter, t, fs_orig):

    fill = np.zeros((len(signal), upsamp_factor - 1))
    signal_up = np.c_[signal, fill].reshape(len(signal) * upsamp_factor,)

    signal_out = fftconvolve(signal_up, fir_lp_filter)

    filter_offset = (len(fir_lp_filter) - 1) / 2

    dt_up = 1.0 / (fs_orig * upsamp_factor)
    t0 = t[0] - filter_offset * dt_up
    t_out = np.arange(0, len(signal_out)) * dt_up + t0

    return signal_out, t_out
