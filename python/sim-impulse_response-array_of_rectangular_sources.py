#!/usr/bin/env python3
# This file is in the public domain.

# Simulate the impulse response at a point due to an array of rectangular
# areas.

import numpy as np
import matplotlib.pyplot as plt
import acoustic_field.excitation as exc
from util import arrayutil
from scipy import signal
import sys

USE_NUMERIC_METHOD = False
if USE_NUMERIC_METHOD:
    import acoustic_field.numeric_rectangular_source_acoustic_field as n_af
else:
    import acoustic_field.analytic_rectangular_source_acoustic_field as a_af

#==============================================================================

USE_FOCUS = True
#---
FOCUS_X = 0.0 # m
FOCUS_Y = 0.0 # m
FOCUS_Z = 40.0e-3 # m
#---
#FOCUS_X = 40.0e-3 * np.sin(np.deg2rad(20.0)) # m
#FOCUS_Y = 0.0 # m
#FOCUS_Z = 40.0e-3 * np.cos(np.deg2rad(20.0)) # m
#---

ARRAY_PITCH_X = 0.6e-3 # m
NUM_ARRAY_ELEM_X = 32
ARRAY_PITCH_Y = 13.0e-3 # m
NUM_ARRAY_ELEM_Y = 1

ELEM_WIDTH  = 0.5e-3 # m
ELEM_HEIGHT = 12.0e-3 # m
C = 1500.0 # m/s
DENSITY = 1000.0 # kg/m3
CENTER_FREQ = 5.0e6 # Hz
MAX_FREQ = CENTER_FREQ * 2.0
nyquist_rate = MAX_FREQ * 2.0
# None: Use default.
EXCITATION_NUM_PERIODS = None

if USE_NUMERIC_METHOD:
    SAMPLE_RATE = nyquist_rate * 8.0 # Hz
    SUB_ELEM_SIZE = C / (nyquist_rate * 4.0) # m
else:
    SAMPLE_RATE = nyquist_rate * 32.0 # Hz
    #MIN_ELEM_EDGE_DIVISOR = 10.0
    MIN_ELEM_EDGE_DIVISOR = None # disable

#==============================================================================

def get_array_impulse_response(x, y, z):

    num_elem = len(elem_x)

    offset_list = []
    h_list = []
    i_min = sys.maxsize
    i_max = 0 # the index after the last
    # Obtain the impulse responses.
    for ie in range(num_elem):
        xe = elem_x[ie]
        ye = elem_y[ie]

        offset, h = af.get_impulse_response(x - xe, y - ye, z)
        offset_with_delay = offset + int(round(focus_delay[ie]))
        offset_list.append(offset_with_delay)
        h_list.append(h)
        i_min = min(i_min, offset_with_delay)
        i_max = max(i_max, offset_with_delay + len(h))
    # Accumulate the impulse responses.
    h_sum = np.zeros(i_max - i_min)
    for ie in range(num_elem):
        i_begin = offset_list[ie]
        i_end = i_begin + len(h_list[ie])
        h_sum[i_begin-i_min:i_end-i_min] += h_list[ie]

    print("x={} y={} z={}".format(x, y, z))
    return i_min, h_sum

#==============================================================================

def plot_h_p(x, y, z):
    offset, h = get_array_impulse_response(x, y, z)
    t_h = arrayutil.get_time_sequence(len(h), SAMPLE_RATE, offset / SAMPLE_RATE)

    half_w = ELEM_WIDTH / 2

    plt.figure()
    plt.plot(C * t_h / half_w, h / C)
    plt.grid(True)
    plt.xlabel("ct/a")
    plt.ylabel("h/c")
    plt.title("h x={} y={} z={}".format(x, y, z))

    p = DENSITY * signal.convolve(h, dvdt)
    t_p = arrayutil.get_time_sequence(len(p), SAMPLE_RATE, offset / SAMPLE_RATE)

    plt.figure()
    plt.plot(C * t_p / half_w, p)
    plt.grid(True)
    plt.xlabel("ct/a")
    plt.ylabel("p")
    plt.title("p x={} y={} z={}".format(x, y, z))

#==============================================================================

def plot_excitation(s):
    plt.figure()
    plt.plot(arrayutil.get_time_sequence(len(s), SAMPLE_RATE), s)
    plt.grid(True)
    plt.title("Excitation waveform")
    plt.xlabel("t (s)")

#==============================================================================

if __name__ == '__main__':

    # Calculate the center of each element.
    end_w = (NUM_ARRAY_ELEM_X - 1) * 0.5 * ARRAY_PITCH_X
    end_h = (NUM_ARRAY_ELEM_Y - 1) * 0.5 * ARRAY_PITCH_Y
    elem_x = np.zeros(NUM_ARRAY_ELEM_X * NUM_ARRAY_ELEM_Y)
    elem_y = np.zeros(NUM_ARRAY_ELEM_X * NUM_ARRAY_ELEM_Y)
    for iy in range(NUM_ARRAY_ELEM_Y):
        for ix in range(NUM_ARRAY_ELEM_X):
            elem_x[iy * NUM_ARRAY_ELEM_X + ix] = ix * ARRAY_PITCH_X - end_w
            elem_y[iy * NUM_ARRAY_ELEM_X + ix] = iy * ARRAY_PITCH_Y - end_h

    # Plot the elements.
    half_w = ELEM_WIDTH / 2
    half_h = ELEM_HEIGHT / 2
    plt.figure()
    for iy in range(NUM_ARRAY_ELEM_Y):
        for ix in range(NUM_ARRAY_ELEM_X):
            x = elem_x[iy * NUM_ARRAY_ELEM_X + ix]
            y = elem_y[iy * NUM_ARRAY_ELEM_X + ix]
            plt.vlines(np.array([x - half_w, x + half_w]), y - half_h, y + half_h)
            plt.hlines(np.array([y - half_h, y + half_h]), x - half_w, x + half_w)
    plt.grid(True)
    plt.title("Array elements")

    # Calculate the focalization delays.
    if USE_FOCUS:
        focus_dt = np.sqrt((FOCUS_X - elem_x)**2 + (FOCUS_Y - elem_y)**2 + FOCUS_Z**2) / C
        max_dt = np.max(focus_dt)
        focus_delay = (max_dt - focus_dt) * SAMPLE_RATE
    else:
        focus_delay = np.zeros(NUM_ARRAY_ELEM_X * NUM_ARRAY_ELEM_Y)

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
    #dvdt = exc.hdf5_waveform("data/ref_pulse-40MHz.h5", "ascan", 40.0e6, SAMPLE_RATE)
    #plot_excitation(dvdt)

    plot_h_p(0.25e-3, 0.0, 40.0e-3)

    plt.show()
