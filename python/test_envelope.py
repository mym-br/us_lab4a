#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert
from util import arrayutil

from util.plotutil import PlotConfig
pcfg = PlotConfig('fig/test_envelope')
pcfg.set_tick_labels_padding(6.0)

FS = 4.0 * 10.0e6
FC = 5.0e6
NUM_PERIODS = 3.0
MARGIN = 20
ENV_PADDING = 16 * 1024
OFFSET = 0.5
TRUNC_LEN = 16

#------------------------------------------------------------------------------

def calc_envelope(s, offset):
    env = np.abs(hilbert(np.r_[s, offset * np.ones((ENV_PADDING,))]))
    env = env[:-ENV_PADDING]
    return env



period = 1.0 / FC

t0 = np.arange(0.0, NUM_PERIODS * period, 1.0 / FS)
s0 = np.sin(2.0 * np.pi * FC * t0)
s1 = -np.cos(2.0 * np.pi * (FC / NUM_PERIODS) * t0) + 1.0

plt.figure()
plt.plot(t0, s0)

plt.figure()
plt.plot(t0, s1)

s2 = s1 * s0
s2 = np.r_[np.zeros((MARGIN,)), s2, np.zeros((MARGIN,))]
print('s2.shape:', s2.shape)
t2 = arrayutil.get_time_sequence(len(s2), FS, 0) * 1e6
s2_env = calc_envelope(s2, 0.0)

plt.figure(figsize=(6, 3))
plt.plot(t2, s2, '--k', label='sinal')
plt.plot(t2, s2_env, 'k', label='envelope')
plt.grid(True)
if pcfg.t3():
    plt.legend(loc='upper right', labelspacing=0.0, fontsize=12, ncol=1)
    plt.xlabel(r'$\mathrm{tempo}~(\si{\micro\second})$')
    plt.ylabel('valor')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/envelope_normal.pdf')

s3 = s2 + OFFSET
s3_env = calc_envelope(s3, OFFSET)

plt.figure(figsize=(6, 3))
plt.plot(t2, s3, '--k', label='sinal')
plt.plot(t2, s3_env, 'k', label='envelope')
plt.grid(True)
if pcfg.t3():
    plt.legend(loc='upper right', labelspacing=0.0, fontsize=12, ncol=1)
    plt.xlabel(r'$\mathrm{tempo}~(\si{\micro\second})$')
    plt.ylabel('valor')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/envelope_dc.pdf')

mask = np.r_[np.zeros((MARGIN + TRUNC_LEN,)), np.ones((len(s0) - TRUNC_LEN + MARGIN,))]
s4 = s2 * mask
s4_env = calc_envelope(s4, 0.0)

plt.figure(figsize=(6, 3))
plt.plot(t2, s4, '--k', label='sinal')
plt.plot(t2, s4_env, 'k', label='envelope')
plt.grid(True)
if pcfg.t3():
    plt.legend(loc='upper right', labelspacing=0.0, fontsize=12, ncol=1)
    plt.xlabel(r'$\mathrm{tempo}~(\si{\micro\second})$')
    plt.ylabel('valor')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/envelope_trunc.pdf')

s5 = np.minimum(s2, 1.0)
s5 = np.maximum(s5, -1.0)
s5_env = calc_envelope(s5, 0.0)

plt.figure()
plt.plot(t2, s5, '--k')
plt.plot(t2, s5_env, 'k')
plt.grid(True)

plt.show()
