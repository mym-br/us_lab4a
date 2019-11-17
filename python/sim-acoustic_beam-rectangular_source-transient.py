#!/usr/bin/env python3
# This file is in the public domain.

# Simulate the acoustic beam of a rectangular area.
# Use multiprocessing.

import numpy as np
import matplotlib.pyplot as plt
import acoustic_field.excitation as exc
from util import arrayutil
from scipy import signal
import time
import multiprocessing as mproc

USE_NUMERIC_METHOD = True
if USE_NUMERIC_METHOD:
    import acoustic_field.numeric_rectangular_source_acoustic_field as n_af
else:
    import acoustic_field.analytic_rectangular_source_acoustic_field as a_af

#==============================================================================

BEAM_DISTANCE = 1.0 # m
BEAM_THETA_Y_STEP = 0.5 # degree
BEAM_THETA_Y_MAX = 90.0 # degree
BEAM_THETA_X_STEP = 0.01 # degree
BEAM_THETA_X_MAX = 10.0 # degree

ELEM_WIDTH  = 0.5e-3 # m
ELEM_HEIGHT = 12.0e-3 # m
C = 1500.0 # m/s
CENTER_FREQ = 5.0e6 # Hz
MAX_FREQ = CENTER_FREQ * 2.0
nyquist_rate = MAX_FREQ * 2.0
# None: Use default.
EXCITATION_NUM_PERIODS = None

if USE_NUMERIC_METHOD:
    SAMPLE_RATE = nyquist_rate * 8.0 # Hz
    SUB_ELEM_SIZE = C / (nyquist_rate * 4.0) # m
else:
    SAMPLE_RATE = nyquist_rate * 512.0 # Hz
    #MIN_ELEM_EDGE_DIVISOR = 10.0
    MIN_ELEM_EDGE_DIVISOR = None # disable

SHOW_POLAR = False

#==============================================================================

def process(acoustic_field, dvdt, len_theta_y,
            grid_x, grid_y, grid_z,
            num_jobs, job_queue, result_queue):
    for job in iter(job_queue.get, "END"):

        beam_theta_x = np.zeros(len_theta_y)
        for iy in range(len_theta_y):
            offset, h = acoustic_field.get_impulse_response(grid_x[job.ix, iy],
                                                            grid_y[job.ix, iy],
                                                            grid_z[job.ix, iy])
            h_dvdt = signal.convolve(h, dvdt)
            beam_theta_x[iy] = np.max(h_dvdt) - np.min(h_dvdt)

        print("ix: {} < {}".format(job.ix, num_jobs))

        result = Result()
        result.ix = job.ix
        result.beam_theta_x = beam_theta_x
        result_queue.put(result)

#==============================================================================

if __name__ == '__main__':

    class Job: pass
    class Result: pass

    start_time = time.time()

    theta_x_list = np.linspace(0.0, BEAM_THETA_X_MAX, int(np.ceil(BEAM_THETA_X_MAX / BEAM_THETA_X_STEP)) + 1)
    theta_y_list = np.linspace(0.0, BEAM_THETA_Y_MAX, int(np.ceil(BEAM_THETA_Y_MAX / BEAM_THETA_Y_STEP)) + 1)

    grid_theta_x, grid_theta_y = np.meshgrid(theta_x_list, theta_y_list, indexing='ij')

    grid_x = BEAM_DISTANCE * np.sin(np.deg2rad(grid_theta_y))
    rx = BEAM_DISTANCE * np.cos(np.deg2rad(grid_theta_y))
    grid_y = rx * np.sin(np.deg2rad(grid_theta_x))
    grid_z = rx * np.cos(np.deg2rad(grid_theta_x))

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

    beam = np.zeros(grid_x.shape)

    job_queue = mproc.Queue()
    result_queue = mproc.Queue()

    # Create jobs.
    num_jobs = len(theta_x_list)
    for ix in range(num_jobs):
        job = Job()
        job.ix = ix
        job_queue.put(job)
    # Create processes.
    proc_list = []
    for i in range(mproc.cpu_count()):
        proc = mproc.Process(target=process, args=(af,
                                                   dvdt,
                                                   len(theta_y_list),
                                                   grid_x,
                                                   grid_y,
                                                   grid_z,
                                                   num_jobs,
                                                   job_queue,
                                                   result_queue))
        proc_list.append(proc)
        proc.start()
    # Get results.
    for i in range(num_jobs):
        result = result_queue.get()
        beam[result.ix, :] = result.beam_theta_x
    # Signal the processes to end their execution.
    for i in range(mproc.cpu_count()):
        job_queue.put("END")
    # Wait for the end of the processes.
    for proc in proc_list:
        # join() must be called after the result queue is emptied.
        proc.join()

    end_time = time.time()
    print("Processing time: {} s".format(end_time - start_time))

    th_y_step = theta_y_list[1] - theta_y_list[0]
    th_x_step = theta_x_list[1] - theta_x_list[0]
    plt.figure()
    #plt.pcolormesh(grid_theta_y - 0.5 * th_y_step, grid_theta_x - 0.5 * th_x_step, beam)
    plt.pcolormesh(grid_theta_y - 0.5 * th_y_step, grid_theta_x - 0.5 * th_x_step,
                   20.0 * np.log10(beam / np.max(np.abs(beam))))
    plt.grid(True)
    #plt.axis("Equal")
    plt.xlabel("theta_y (degree)")
    plt.ylabel("theta_x (degree)")
    plt.colorbar()

    fig = plt.figure()
    ax1 = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=SHOW_POLAR)
    ax1.set_ylim(-60.0, 0.0)
    ax1.set_yticks(np.arange(-60.0, 0.1, 10.0))
    if SHOW_POLAR:
        ax1.plot(np.deg2rad(theta_y_list), 20.0 * np.log10(beam[0, :] / np.max(np.abs(beam[0, :]))))
    else:
        ax1.plot(theta_y_list, 20.0 * np.log10(beam[0, :] / np.max(np.abs(beam[0, :]))))
        plt.xlabel("degree")
        plt.ylabel("dB")
    plt.grid(True)
    plt.title("theta_y")

    fig = plt.figure()
    ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=SHOW_POLAR)
    ax2.set_ylim(-60.0, 0.0)
    ax2.set_yticks(np.arange(-60.0, 0.1, 10.0))
    if SHOW_POLAR:
        ax2.plot(np.deg2rad(theta_x_list), 20.0 * np.log10(beam[:, 0] / np.max(np.abs(beam[:, 0]))))
    else:
        ax2.plot(theta_x_list, 20.0 * np.log10(beam[:, 0] / np.max(np.abs(beam[:, 0]))))
        plt.xlabel("degree")
        plt.ylabel("dB")
    plt.grid(True)
    plt.title("theta_x")

    plt.show()
