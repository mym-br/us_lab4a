#!/usr/bin/env python3
# This file is in the public domain.

# Simulate the acoustic field of an array of rectangular areas.
# Use multiprocessing.

import numpy as np
import matplotlib.pyplot as plt
import acoustic_field.excitation as exc
from util import arrayutil, hdf5util
from scipy import signal
import time
import sys
import multiprocessing as mproc

USE_NUMERIC_METHOD = False
if USE_NUMERIC_METHOD:
    import acoustic_field.numeric_rectangular_source_acoustic_field as n_af
else:
    import acoustic_field.analytic_rectangular_source_acoustic_field as a_af

_EPS = np.spacing(1)

#==============================================================================

USE_FOCUS = True
#---
#FOCUS_X = 0.0 # m
#FOCUS_Y = 0.0 # m
#FOCUS_Z = 40.0e-3 # m
#---
FOCUS_X = 40.0e-3 * np.sin(np.deg2rad(20.0)) # m
FOCUS_Y = 0.0 # m
FOCUS_Z = 40.0e-3 * np.cos(np.deg2rad(20.0)) # m
#---

ARRAY_PITCH_X = 0.6e-3 # m
NUM_ARRAY_ELEM_X = 32
ARRAY_PITCH_Y = 13.0e-3 # m
NUM_ARRAY_ELEM_Y = 1

ELEM_WIDTH  = 0.5e-3 # m
ELEM_HEIGHT = 12.0e-3 # m
C = 1500.0 # m/s
CENTER_FREQ = 5.0e6 # Hz
MAX_FREQ = CENTER_FREQ * 2.0
nyquist_rate = MAX_FREQ * 2.0
# None: Use default.
EXCITATION_NUM_PERIODS = None

nyq_lambda = C / nyquist_rate
X_MIN = -30.0e-3 # m
X_MAX = 30.0e-3 # m
X_DIV = 1.0
x_step = nyq_lambda / X_DIV
Z_MIN = 1.0e-3 # m
Z_MAX = 60.0e-3 # m
Z_DIV = 1.0
z_step = nyq_lambda / Z_DIV

Y = 0.0 # m

if USE_NUMERIC_METHOD:
    SAMPLE_RATE = nyquist_rate * 8.0 # Hz
    SUB_ELEM_SIZE = C / (nyquist_rate * 2.0) # m
else:
    SAMPLE_RATE = nyquist_rate * 32.0 # Hz
    #MIN_ELEM_EDGE_DIVISOR = 10.0
    MIN_ELEM_EDGE_DIVISOR = None # disable

#==============================================================================

def process(acoustic_field, dvdt, z_list, elem_x, elem_y, focus_delay, num_jobs, job_queue, result_queue):
    for job in iter(job_queue.get, "END"):

        num_elem = len(elem_x)

        image_line = np.zeros(len(z_list))
        for iz, z in enumerate(z_list):
            offset_list = []
            h_list = []
            i_min = sys.maxsize
            i_max = 0 # the index after the last
            # Obtain the impulse responses.
            for ie in range(num_elem):
                xe = elem_x[ie]
                ye = elem_y[ie]

                offset, h = acoustic_field.get_impulse_response(job.x - xe, job.y - ye, z)
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

            h_dvdt = signal.convolve(h_sum, dvdt)
            image_line[iz] = np.max(h_dvdt) - np.min(h_dvdt)

        print("ix: {} < {}".format(job.ix, num_jobs))

        result = Result()
        result.ix = job.ix
        result.image_line = image_line
        result_queue.put(result)

#==============================================================================

if __name__ == '__main__':

    class Job: pass
    class Result: pass

    start_time = time.time()

    x_list = np.arange(X_MIN, X_MAX + _EPS, x_step)
    z_list = np.arange(Z_MIN, Z_MAX + _EPS, z_step)

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
    #v = exc.hdf5_waveform("data/ref_pulse-40MHz.h5", "ascan", 40.0e6, SAMPLE_RATE)
    plt.figure()
    plt.plot(arrayutil.get_time_sequence(len(v), SAMPLE_RATE), v)
    plt.grid(True)
    plt.title("Excitation waveform")
    plt.xlabel("t (s)")

    dvdt = arrayutil.central_diff(v, SAMPLE_RATE)

    grid_x, grid_z = np.meshgrid(x_list, z_list, indexing='ij')
    image = np.zeros(grid_x.shape)

    job_queue = mproc.Queue()
    result_queue = mproc.Queue()

    # Create jobs.
    num_jobs = len(x_list)
    for ix in range(num_jobs):
        job = Job()
        job.ix = ix
        job.y = Y
        job.x = x_list[ix]
        job_queue.put(job)
    # Create processes.
    proc_list = []
    for i in range(mproc.cpu_count()):
        proc = mproc.Process(target=process, args=(af,
                                                   dvdt,
                                                   z_list,
                                                   elem_x,
                                                   elem_y,
                                                   focus_delay,
                                                   num_jobs,
                                                   job_queue,
                                                   result_queue))
        proc_list.append(proc)
        proc.start()
    # Get results.
    for i in range(num_jobs):
        result = result_queue.get()
        image[result.ix, :] = result.image_line
    # Signal the processes to end their execution.
    for i in range(mproc.cpu_count()):
        job_queue.put("END")
    # Wait for the end of the processes.
    for proc in proc_list:
        # join() must be called after the result queue is emptied.
        proc.join()

    end_time = time.time()
    print("Processing time: {} s".format(end_time - start_time))

    plt.figure()
    plt.pcolormesh(grid_z - 0.5 * z_step, grid_x - 0.5 * x_step, image)
    plt.grid(True)
    #plt.axis("Equal")
    plt.xlabel("z (m)")
    plt.ylabel("x (m)")

    hdf5util.write_ndarray(image, "output/image_value.h5", "value")
    hdf5util.write_ndarray(grid_x, "output/image_x.h5", "x")
    hdf5util.write_ndarray(grid_z, "output/image_z.h5", "z")

    plt.show()
