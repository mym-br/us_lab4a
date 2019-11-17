#!/usr/bin/env python3
# This file is in the public domain.

# Simulate the impulse response at a point due to a rectangular area.

import matplotlib.pyplot as plt
import acoustic_field.excitation as exc
from util import arrayutil
from scipy import signal

USE_NUMERIC_METHOD = True
if USE_NUMERIC_METHOD:
    import acoustic_field.numeric_rectangular_source_acoustic_field as n_af
else:
    import acoustic_field.analytic_rectangular_source_acoustic_field as a_af

#==============================================================================

ELEM_WIDTH  = 10.0e-3 # m
ELEM_HEIGHT = 16.0e-3 # m
C = 1500.0 # m/s
DENSITY = 1000.0 # kg/m3
CENTER_FREQ = 5.0e6 # Hz
MAX_FREQ = CENTER_FREQ * 2.0
nyquist_rate = MAX_FREQ * 2.0
# None: Use default.
EXCITATION_NUM_PERIODS = None

if USE_NUMERIC_METHOD:
    SAMPLE_RATE = nyquist_rate * 4.0 # Hz
    SUB_ELEM_SIZE = C / (nyquist_rate * 4.0) # m
else:
    SAMPLE_RATE = nyquist_rate * 8.0 # Hz
    MIN_ELEM_EDGE_DIVISOR = 10.0
    #MIN_ELEM_EDGE_DIVISOR = None # disable

#==============================================================================

def plot_h_p(x, y, z):
    offset, h = af.get_impulse_response(x, y, z)
    t_h = arrayutil.get_time_sequence(len(h), SAMPLE_RATE, offset / SAMPLE_RATE)

    half_w = ELEM_WIDTH / 2

    plt.figure()
    plt.plot(C * t_h / half_w, h / C)
    plt.grid(True)
    plt.xlabel("ct/a")
    plt.ylabel("h/c")
    plt.title("h x/a={} y/a={} z/a={}".format(x / half_w, y / half_w, z / half_w))

    p = DENSITY * signal.convolve(h, dvdt)
    t_p = arrayutil.get_time_sequence(len(p), SAMPLE_RATE, offset / SAMPLE_RATE)

    plt.figure()
    plt.plot(C * t_p / half_w, p)
    plt.grid(True)
    plt.xlabel("ct/a")
    plt.ylabel("p")
    plt.title("p x/a={} y/a={} z/a={}".format(x / half_w, y / half_w, z / half_w))

#==============================================================================

def plot_excitation(s):
    plt.figure()
    plt.plot(arrayutil.get_time_sequence(len(s), SAMPLE_RATE), s)
    plt.grid(True)
    plt.title("Excitation waveform")
    plt.xlabel("t (s)")

#==============================================================================

if __name__ == '__main__':

    if USE_NUMERIC_METHOD:
        af = n_af.NumericRectangularSourceAcousticField(ELEM_WIDTH, ELEM_HEIGHT, SAMPLE_RATE, C, SUB_ELEM_SIZE)
        af.plot_sub_elements()
    else:
        af = a_af.AnalyticRectangularSourceAcousticField(ELEM_WIDTH, ELEM_HEIGHT, SAMPLE_RATE, C, MIN_ELEM_EDGE_DIVISOR)

    #v = exc.waveform1(CENTER_FREQ, SAMPLE_RATE, EXCITATION_NUM_PERIODS)
    v = exc.waveform2a(CENTER_FREQ, SAMPLE_RATE, EXCITATION_NUM_PERIODS)
    #v = exc.waveform2b(CENTER_FREQ, SAMPLE_RATE, EXCITATION_NUM_PERIODS)
    plot_excitation(v)
    dvdt = arrayutil.central_diff(v, SAMPLE_RATE)

    half_w = ELEM_WIDTH / 2

    plot_h_p(2.0 * half_w, 2.0 * half_w, 5.0 * half_w)
    plot_h_p(0.4 * half_w, 2.0 * half_w, 5.0 * half_w)
    plot_h_p(2.0 * half_w, 0.4 * half_w, 5.0 * half_w)
    plot_h_p(0.4 * half_w, 0.4 * half_w, 5.0 * half_w)

    plt.show()
