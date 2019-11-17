#!/usr/bin/env python3
# This file is in the public domain.

import numpy as np

def _fade_in(n):
    n = float(n)
    return np.arange(n) / n

def _fade_out(n):
    return _fade_in(n)[::-1]

def get_rectangular_window(n):
    return np.ones((n,))

def get_trapezoidal_window(n_fade_in, n_top, n_fade_out):
    return np.r_[_fade_in(n_fade_in), np.ones((n_top,)), _fade_out(n_fade_out)]

# The first and the last zeroed elements are not included in the returned array.
def get_trapezoidal2_window(n_fade_in, n_top, n_fade_out):
    return np.r_[_fade_in(n_fade_in + 1)[1:], np.ones((n_top,)), _fade_out(n_fade_out + 1)[:-1]]

def get_triangular_window(n):
    n2 = 0.5 * (n - 1)
    if n % 2 == 0:
        return np.r_[_fade_in(n2), _fade_out(n2)]
    else:
        return np.r_[_fade_in(n2), 1.0, _fade_out(n2)]

# The first and the last zeroed elements are not included in the returned array.
def get_triangular2_window(n):
    n2 = 0.5 * (n - 1)
    if n % 2 == 0:
        return np.r_[_fade_in(n2 + 1)[1:], _fade_out(n2 + 1)[:-1]]
    else:
        return np.r_[_fade_in(n2 + 1)[1:], 1.0, _fade_out(n2 + 1)[:-1]]

def get_hanning_window(n):
    return 0.5 - 0.5 * np.cos(2.0 * np.pi * np.linspace(0.0, 1.0, n))

# The first and the last zeroed elements are not included in the returned array.
def get_hanning2_window(n):
    return 0.5 - 0.5 * np.cos(2.0 * np.pi * np.linspace(0.0, 1.0, n + 2)[1:-1])

def get_hamming_window(n):
    return 0.54 - 0.46 * np.cos(2.0 * np.pi * np.linspace(0.0, 1.0, n))

def get_blackman_window(n):
    m = np.linspace(0.0, 1.0, n)
    return 0.42 - 0.5 * np.cos(2.0 * np.pi * m) + 0.08 * np.cos(4.0 * np.pi * m)

# The first and the last zeroed elements are not included in the returned array.
def get_blackman2_window(n):
    m = np.linspace(0.0, 1.0, n + 2)[1:-1]
    return 0.42 - 0.5 * np.cos(2.0 * np.pi * m) + 0.08 * np.cos(4.0 * np.pi * m)



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    plt.stem(get_blackman2_window(15))
    plt.grid(True)
    plt.show()
