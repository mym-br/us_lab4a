#!/usr/bin/env python3
# This file is in the public domain.

# Simulate the acoustic field of a rectangular area.
# Use multiprocessing.

import numpy as np
import matplotlib.pyplot as plt
import acoustic_field.excitation as exc
from util import arrayutil
from scipy import signal
import time
import multiprocessing as mproc

USE_NUMERIC_METHOD = False
if USE_NUMERIC_METHOD:
    import acoustic_field.numeric_rectangular_source_acoustic_field as n_af
else:
    import acoustic_field.analytic_rectangular_source_acoustic_field as a_af

_EPS = np.spacing(1)

#==============================================================================

ELEM_WIDTH  = 10.0e-3 # m
ELEM_HEIGHT = 16.0e-3 # m
C = 1500.0 # m/s
CENTER_FREQ = 5.0e6 # Hz
MAX_FREQ = CENTER_FREQ * 2.0
nyquist_rate = MAX_FREQ * 2.0
# None: Use default.
EXCITATION_NUM_PERIODS = None

nyq_lambda = C / nyquist_rate
X_MIN = -20.0e-3 # m
X_MAX = 20.0e-3 # m
X_DIV = 2.0
x_step = nyq_lambda / X_DIV
Z_MIN = 1.0e-3 # m
Z_MAX = 60.0e-3 # m
Z_DIV = 2.0
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

def process(acoustic_field, dvdt, z_list, num_jobs, job_queue, result_queue):
    for job in iter(job_queue.get, "END"):

        image_line = np.zeros(len(z_list))
        for iz, z in enumerate(z_list):
            offset, h = acoustic_field.get_impulse_response(job.x, job.y, z)
            h_dvdt = signal.convolve(h, dvdt)
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

    if USE_NUMERIC_METHOD:
        af = n_af.NumericRectangularSourceAcousticField(ELEM_WIDTH, ELEM_HEIGHT, SAMPLE_RATE, C, SUB_ELEM_SIZE)
        af.plot_sub_elements()
    else:
        af = a_af.AnalyticRectangularSourceAcousticField(ELEM_WIDTH, ELEM_HEIGHT, SAMPLE_RATE, C, MIN_ELEM_EDGE_DIVISOR)

    #v = exc.waveform1(CENTER_FREQ, SAMPLE_RATE, EXCITATION_NUM_PERIODS)
    v = exc.waveform2a(CENTER_FREQ, SAMPLE_RATE, EXCITATION_NUM_PERIODS)
    #v = exc.waveform2b(CENTER_FREQ, SAMPLE_RATE, EXCITATION_NUM_PERIODS)
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

    plt.show()
