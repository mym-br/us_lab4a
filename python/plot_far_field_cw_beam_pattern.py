#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np

from util.plotutil import PlotConfig
pcfg = PlotConfig('fig/plot_far_field_cw_beam_pattern')
pcfg.set_tick_labels_padding(6.0)

MIN_ANGLE = -90.0
MAX_ANGLE = 90.0
NUM_ANGLES = 10000

FC = 5.0e6
PITCH = 0.15e-3
C = 1500.0
NUM_ELEM = 32

MIN_LOG = -40.0

#------------------------------------------------------------------------------

lambda_c = C / FC
print('PITCH / lambda_c:', PITCH / lambda_c)

angle_list_deg = np.linspace(MIN_ANGLE, MAX_ANGLE, NUM_ANGLES)
angle_list = np.deg2rad(angle_list_deg)

fe = -np.sin(angle_list) / lambda_c
omega_n = 2.0 * np.pi * fe * PITCH

s = np.zeros(angle_list.shape, dtype=complex)
for m in np.arange(0, NUM_ELEM, dtype=float):
    s += np.exp(-1j * m * omega_n)

s2 = np.abs(s)
s2 /= s2.max()
s2 = np.maximum(10.0**(MIN_LOG / 20.0), s2)

plt.figure(figsize=(6, 4))
plt.plot(angle_list_deg, 20.0 * np.log10(s2), 'k')
plt.xlim(-90.0, 90.0)
plt.grid(True)
if pcfg.t3():
    plt.xlabel(r'$\theta_\mathrm{f} (\si{\degree})$')
    plt.ylabel('amplitude (dB)')
    plt.autoscale(enable=True, axis='x', tight=True)
    plt.tight_layout(pad=0.2) # padding between the x and y labels and the figure borders
    plt.savefig(pcfg.plot_dir + '/beam_pattern.pdf')

plt.show()
