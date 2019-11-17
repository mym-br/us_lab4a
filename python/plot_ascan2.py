#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from scipy.fftpack import fft
from scipy.signal import hilbert
from util import arrayutil, hdf5util, interpolation

from util.plotutil import PlotConfig
pcfg = PlotConfig('fig/plot_ascan2')
pcfg.set_tick_labels_padding(6.0)

ASCAN_FILE_NAME = '../data/ref_pulse-40MHz.h5'
ASCAN_DATASET_NAME = 'ascan'
FS = 40.0e6

UPSAMPLING_FACTOR = 16
UPSAMP_FILTER_TOLERANCE = 0.0001
UPSAMP_FILTER_TRANSITION_WIDTH = 0.2 # fraction of the original bandwidth

MAX_FREQ = 20e6
ENV_PADDING = 1024
FREQ_RESP_THRESHOLD = 0.1
FREQ_RESP_PADDING = 32 * 1024

#------------------------------------------------------------------------------

fs = FS * UPSAMPLING_FACTOR

ascan = hdf5util.read_to_ndarray(ASCAN_FILE_NAME, ASCAN_DATASET_NAME)
t = arrayutil.get_time_sequence(len(ascan), FS, 0)

plt.figure(figsize=(6, 3))
plt.plot(t, ascan / np.abs(ascan).max(), 'k')
plt.grid(True)
if pcfg.t1():
    plt.title('A-scan')
    plt.xlabel('time (s)')
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
elif pcfg.t2():
    plt.xlabel(r'$\mathrm{time} (\si{\second})$')
    plt.ylabel('normalized value')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/ascan.pdf')
elif pcfg.t3():
    plt.xlabel(r'$\mathrm{tempo} (\si{\second})$')
    plt.ylabel('valor normalizado')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/ascan.pdf')

lp_filter = interpolation.upsampling_filter(UPSAMPLING_FACTOR,
                                            UPSAMP_FILTER_TRANSITION_WIDTH,
                                            UPSAMP_FILTER_TOLERANCE,
                                            plot=False)
ascan_r, t_r = interpolation.interpolate(ascan,
                                         UPSAMPLING_FACTOR,
                                         lp_filter,
                                         t,
                                         FS)
env = np.abs(hilbert(np.r_[ascan_r, np.zeros((ENV_PADDING,))]))
env = env[:-ENV_PADDING]
t_env = arrayutil.get_time_sequence(len(env), fs, t_r[0])

plt.figure(figsize=(6, 3))
plt.plot(t_r, ascan_r / env.max(), 'k', label='A-scan')
plt.plot(t_env, env / env.max(), 'k--', label='Envelope')
plt.legend(loc='upper right', fontsize=12, labelspacing=0.0)
plt.grid(True)
plt.ylabel('normalized value')
if pcfg.t1():
    plt.title('A-scan (interpolated)')
    plt.xlabel('time (s)')
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
elif pcfg.t2():
    plt.xlabel(r'$\mathrm{time} (\si{\second})$')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/ascan_r.pdf')

freq_resp_padding = np.zeros((FREQ_RESP_PADDING,))
freq_resp = np.abs(fft(np.r_[ascan_r, freq_resp_padding]))
freq_resp /= freq_resp.max()
freq = np.arange(len(freq_resp)) * fs / len(freq_resp)
max_idx = int((MAX_FREQ / fs) * len(freq))

plt.figure(figsize=(6, 3))
plt.plot(freq[:max_idx] * 1e-6, freq_resp[:max_idx], 'k')
plt.grid(True)
plt.ylabel('normalized value')
if pcfg.t1():
    plt.title('Freq. response of ref. pulse')
    plt.xlabel('f (MHz)')
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
elif pcfg.t2():
    plt.xlabel(r'$\mathrm{frequency} (\si{\mega\hertz})$')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/freq_resp.pdf')

freq_resp_env = np.abs(fft(np.r_[env, freq_resp_padding]))
freq_resp_env /= freq_resp_env.max()
freq_env = np.arange(len(freq_resp_env)) * fs / len(freq_resp_env)
max_idx_env = int((MAX_FREQ / fs) * len(freq_env))

plt.figure(figsize=(6, 3))
plt.plot(freq_env[:max_idx_env] * 1e-6, freq_resp_env[:max_idx_env], 'k')
plt.grid(True)
plt.ylabel('normalized value')
if pcfg.t1():
    plt.title('Freq. response of ref. pulse (envelope)')
    plt.xlabel('f (MHz)')
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
elif pcfg.t2():
    plt.xlabel(r'$\mathrm{frequency} (\si{\mega\hertz})$')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/freq_resp_env.pdf')

min_freq = freq[freq_resp[:max_idx] >= FREQ_RESP_THRESHOLD][0]
max_freq = freq[freq_resp[:max_idx] >= FREQ_RESP_THRESHOLD][-1]
print('min_freq:', min_freq, 'max_freq:', max_freq)
min_freq_env = freq_env[freq_resp_env[:max_idx_env] >= FREQ_RESP_THRESHOLD][0]
max_freq_env = freq_env[freq_resp_env[:max_idx_env] >= FREQ_RESP_THRESHOLD][-1]
print('min_freq_env:', min_freq_env, 'max_freq_env:', max_freq_env)

plt.show()
