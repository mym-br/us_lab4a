#!/usr/bin/env python3
# This file is in the public domain.

import matplotlib.pyplot as plt
import numpy as np
from util import hdf5util
from scipy.signal import kaiser_beta, kaiserord



# tolerance_db > 0.0
def beta(tolerance_db):
    if tolerance_db > 50.0:
        beta = 0.1102 * (tolerance_db - 8.7)
    elif tolerance_db > 21.0:
        k = tolerance_db - 21.0
        beta = 0.5842 * k**0.4 + 0.07886 * k
    else:
        beta = 0.0
    return beta

# tolerance_db > 0.0
# transition_width: 1.0 --> fs/2
def size(tolerance_db, transition_width):
    # Scipy and Matlab use 7.95.
    # Oppenheim and Octave use 8.0.
    m = np.ceil((tolerance_db - 7.95) / (2.285 * transition_width * np.pi) + 1.0)
    return int(m)

# 0 < tolerance <= 0.3
# transition_width: 1.0 --> fs/2
def param(tolerance, transition_width):
    tol_db = -20.0 * np.log10(tolerance)
    window_size = size(tol_db, transition_width)
    window_beta = beta(tol_db)
    return window_size, window_beta



tol_list = np.arange(1e-6, 0.2, 1e-4)
tol_db_list = -20.0 * np.log10(tol_list)
beta_list = np.empty_like(tol_db_list)
for i, tol in enumerate(np.nditer(tol_db_list)):
    beta_list[i] = beta(tol)

beta_list_ref = np.empty_like(tol_db_list)
for i, tol in enumerate(np.nditer(tol_db_list)):
    beta_list_ref[i] = kaiser_beta(tol)

print('max abs beta error:', np.abs(beta_list - beta_list_ref).max())
hdf5util.write_ndarray(data=tol_db_list, file_path='kaiser_tol_db.h5', dataset_name='v')
hdf5util.write_ndarray(data=beta_list_ref, file_path='kaiser_beta.h5', dataset_name='v')

plt.figure()
plt.plot(tol_db_list, beta_list, label='local')
plt.plot(tol_db_list, beta_list_ref, label='scipy')
plt.title('Beta list')
plt.legend(loc='upper right', labelspacing=0.2)



trans_width_list = np.array([0.05, 0.1, 0.5])
size_matrix = np.zeros((len(trans_width_list), len(tol_list)), dtype=int)
size_matrix_ref = np.empty_like(size_matrix)
for i, trans_width in enumerate(np.nditer(trans_width_list)):
    plt.figure()
    for j, tol in enumerate(np.nditer(tol_list)):
        w_size, w_beta = param(tol, trans_width)
        size_matrix[i, j] = w_size

        w_size_ref, w_beta_ref = kaiserord(-20.0 * np.log10(tol), trans_width)
        size_matrix_ref[i, j] = w_size_ref
    plt.plot(tol_db_list, size_matrix[i, :], label='local')
    plt.plot(tol_db_list, size_matrix_ref[i, :], label='scipy')
    plt.title('Window size (transition width: ' + str(trans_width))
    plt.legend(loc='upper right', labelspacing=0.2)

    print('max abs size error (trans. width = ', trans_width, '):', np.abs(size_matrix[i, :] - size_matrix_ref[i, :]).max())

    hdf5util.write_ndarray(data=size_matrix_ref[i, :].astype(float),
                           file_path='kaiser_size-trans_width_' + str(trans_width) + '.h5',
                           dataset_name='v')

plt.show()
