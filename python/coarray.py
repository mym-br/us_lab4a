#!/usr/bin/env python3
# This file is in the public domain.

# Calculate the coarray for a linear array.
# Plot the directivity (CW, far-field).

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import convolve
import math
import matplotlib as mpl
mpl.style.use('classic')

import util.windowfunction as wf

#==============================================================================

# If True:  The elements are represented by rectangles infinitely long in
#           the y-axis.
# If False: The elements are represented by lines.
USE_RECT_ELEM = True

F = 5.0e6 # Hz
C = 1500.0 # m/s
PITCH = 0.6e-3 # m
MIN_DECIBEL = -100.0
MIN_ANGLE = -90.0 # degree
MAX_ANGLE = 90.0 # degree
ANGLE_STEP = 0.01 # degree

if USE_RECT_ELEM:
    # element width = PITCH * ELEM_WIDTH_DIV / (ELEM_WIDTH_DIV + ELEM_SPACING_DIV)
    ELEM_WIDTH_DIV = 50
    # space between elements = PITCH * ELEM_SPACING_DIV / (ELEM_WIDTH_DIV + ELEM_SPACING_DIV)
    ELEM_SPACING_DIV = 10

#==============================================================================

_lambda = C / F
print('lambda:', _lambda)

if USE_RECT_ELEM:
    dx = PITCH / (ELEM_WIDTH_DIV + ELEM_SPACING_DIV)
else:
    dx = PITCH

class AcqStep: pass

def plot_coarray():
    #plt.figure(figsize=(6, 4))
    plt.figure()
    plt.stem(np.arange(len(coarray)), coarray, markerfmt=".")
    plt.title('Coarray')
    plt.xlabel('element')
    plt.ylabel('value')
    plt.grid(True)
    y_range = coarray.max() - coarray.min()
    plt.axis([-5, len(coarray) + 5, -0.1 * y_range, coarray.max() + 0.1 * y_range])
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders

def plot_directivity():
    #plt.figure(figsize=(6, 4))
    plt.figure()

    num_angles = math.ceil((MAX_ANGLE - MIN_ANGLE) / ANGLE_STEP) + 1
    print('num_angles:', num_angles)
    angle_list = np.linspace(MIN_ANGLE, MAX_ANGLE, num_angles)

    d = np.zeros(len(angle_list))
    for i, angle in enumerate(angle_list):
        sin_angle = np.sin(np.deg2rad(angle))
        d[i] = abs(np.sum(coarray * np.exp(1j * (2.0 * np.pi / _lambda) *
                                           np.arange(len(coarray)) * dx *
                                           sin_angle)))
    min_level = 10.0**(MIN_DECIBEL / 20.0)
    d /= d.max()
    d = np.maximum(min_level, d)
    print('len(coarray):', len(coarray))

    plt.plot(angle_list, 20.0 * np.log10(d), '-k')
    plt.grid(True)
    plt.title('Directivity')
    plt.xlabel('theta (degree)')
    plt.ylabel('normalized value (dB)')
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.xticks(ticks=np.arange(-90.0, 91.0, 10.0))

def get_coarray(num_elem, acq_seq, tx_apod, rx_apod):
    if USE_RECT_ELEM:
        array_step = ELEM_WIDTH_DIV + ELEM_SPACING_DIV
        array_size = array_step * (num_elem - 1) + 1
        elem_array = np.ones((ELEM_WIDTH_DIV,))
        coarray = np.zeros(((array_size + ELEM_WIDTH_DIV - 1) * 2 - 1,))
        for i, acq_step in enumerate(acq_seq):
            tx_active = np.zeros(array_size)
            for local_elem in acq_step.tx_elem:
                elem = acq_step.base_elem + local_elem
                if elem < 0: raise ValueError("Negative element index: {}.".format(elem))
                tx_active[elem * array_step] = tx_apod[local_elem]
            rx_active = np.zeros(array_size)
            for local_elem in acq_step.rx_elem:
                elem = acq_step.base_elem + local_elem
                if elem < 0: raise ValueError("Negative element index: {}.".format(elem))
                rx_active[elem * array_step] = rx_apod[local_elem]
            tx_array = convolve(tx_active, elem_array)
            rx_array = convolve(rx_active, elem_array)
            c = convolve(tx_array, rx_array)
            coarray += c
        return coarray
    else:
        coarray = np.zeros((num_elem * 2 - 1,))
        for i, acq_step in enumerate(acq_seq):
            tx_active = np.zeros(num_elem)
            for local_elem in acq_step.tx_elem:
                elem = acq_step.base_elem + local_elem
                if elem < 0: raise ValueError("Negative element index: {}.".format(elem))
                tx_active[elem] = tx_apod[local_elem]
            rx_active = np.zeros(num_elem)
            for local_elem in acq_step.rx_elem:
                elem = acq_step.base_elem + local_elem
                if elem < 0: raise ValueError("Negative element index: {}.".format(elem))
                rx_active[elem] = rx_apod[local_elem]
            c = convolve(tx_active, rx_active)
            coarray += c
        return coarray

def get_coarray_one_way():
    num_elem = 32
    if USE_RECT_ELEM:
        array_step = ELEM_WIDTH_DIV + ELEM_SPACING_DIV
        array_size = array_step * num_elem
        elem_array = np.ones((ELEM_WIDTH_DIV,))
        active = np.zeros(array_size - array_step + 1)
        for elem in range(num_elem):
            active[elem * array_step] = 1.0
        coarray = convolve(active, elem_array)
        return coarray
    else:
        return np.ones(num_elem)

def get_coarray_saft():
    num_elem = 32
    group_size = num_elem
    acq_seq = []
    for i in range(num_elem):
        step = AcqStep()
        step.base_elem = 0
        step.tx_elem = [i]
        step.rx_elem = [i]
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_cyl():
    num_elem = 32
    group_size = num_elem
    acq_seq = []
    step = AcqStep()
    step.base_elem = 0
    step.tx_elem = [15]
    step.rx_elem = np.arange(group_size)
    acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_sta():
    num_elem = 32
    group_size = num_elem
    acq_seq = []
    for i in range(num_elem):
        step = AcqStep()
        step.base_elem = 0
        step.tx_elem = [i]
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_cyl_wave_config13():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 8):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = [15]
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_trapezoidal2_window(15, 2, 15)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_cyl_wave_config9():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 12):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = [15]
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_trapezoidal2_window(8, 16, 8)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_cyl_wave_config7():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 16):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = [15]
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_cyl_wave_config5():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 24):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = [15]
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_sta_config9():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 12):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = np.arange(group_size)
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_trapezoidal2_window(8, 16, 8)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_sta_config7():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 16):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = np.arange(group_size)
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)

def get_coarray_mx_group_sta_config5():
    num_elem = 128
    group_size = 32
    acq_seq = []
    for base_elem in range(0, 97, 24):
        step = AcqStep()
        step.base_elem = base_elem
        step.tx_elem = np.arange(group_size)
        step.rx_elem = np.arange(group_size)
        acq_seq.append(step)
    tx_apod = wf.get_rectangular_window(group_size)
    rx_apod = wf.get_rectangular_window(group_size)
    return get_coarray(num_elem, acq_seq, tx_apod, rx_apod)



#coarray = get_coarray_one_way()
#coarray = get_coarray_saft()
#coarray = get_coarray_cyl()
coarray = get_coarray_sta()
#coarray = get_coarray_mx_group_cyl_wave_config13()
#coarray = get_coarray_mx_group_cyl_wave_config9()
#coarray = get_coarray_mx_group_cyl_wave_config7()
#coarray = get_coarray_mx_group_cyl_wave_config5()
#coarray = get_coarray_mx_group_sta_config9()
#coarray = get_coarray_mx_group_sta_config7()
#coarray = get_coarray_mx_group_sta_config5()

plot_coarray()
plot_directivity()

plt.show()
