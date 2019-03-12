#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from util import hdf5util
from scipy.special import iv

x = np.arange(0.0, 50.0, 0.1)

bessel_i0 = iv(0.0, x)

plt.figure()
plt.plot(x, bessel_i0)
plt.title('bessel I0')

hdf5util.write_ndarray(data=x, file_path='bessel_i0_x.h5', dataset_name='v')
hdf5util.write_ndarray(data=bessel_i0, file_path='bessel_i0.h5', dataset_name='v')

plt.show()
