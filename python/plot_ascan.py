#!/usr/bin/env python3
# This file is in the public domain.

import h5py
import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) != 3:
	sys.exit("Usage: " + sys.argv[0] + " ascan_file_path sampling_freq")

ascan_file_path = sys.argv[1]
fs = float(sys.argv[2])

ascan_dataset_name = 'ascan'
dt = 1.0 / fs

print('Reading A-scan...')
ascan_dataset = h5py.File(ascan_file_path)[ascan_dataset_name]
ascan = np.zeros(ascan_dataset.shape)
ascan_dataset.read_direct(ascan)

t = np.arange(0, ascan.shape[0]) * (1.0 / fs)

print('plot...')
plt.figure(num=1, figsize=(12, 6), dpi=85)
plt.plot(t, ascan)
plt.title('A-scan')
plt.xlabel('t (s)')
plt.ylabel('value')
plt.grid(True)

print('savefig...')
plt.savefig('ascan.png', dpi=300)

print('show...')
plt.show()
