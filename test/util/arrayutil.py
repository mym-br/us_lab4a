import numpy as np

def normalize(a):
    return a * (1.0 / np.abs(a).max())

def get_time_sequence(a, fs, offset):
    return np.arange(len(a)) * (1.0 / fs) + offset

def central_diff_with_time(t, a, fs):
    ts = 1.0 / fs
    a_d = (np.r_[a[:], np.zeros((2,))] - np.r_[np.zeros((2,)), a[:]]) * (0.5 * fs)
    t_d = np.arange(len(a_d)) + (t[0] - ts)
    return t_d, a_d

def central_diff(a, fs):
    a_d = (np.r_[a[:], np.zeros((2,))] - np.r_[np.zeros((2,)), a[:]]) * (0.5 * fs)
    return a_d

def get_dft_freq(n, fs):
    n = float(n)
    return np.arange(n) * (fs / n)

def get_next_power_of_two(n):
    f = 2.0
    n = float(abs(n))
    while n > f:
        f *= 2.0
    return f

# n: size after padding
def pad_with_zeros(a, n):
    if n < len(a):
        raise ValueError('n < len(a)')
    return np.r_[a, np.zeros((n - len(a),))]

def get_cubic_fade_curve(v0, v1, t):
    curve = (2 * t**3 - 3 * t**2 + 1) * v0 + (-2 * t**3 + 3 * t**2) * v1
    return curve

def get_next_fast_even_fft_size(n):
    """Gets the next fast even fft size greater than or equal to n.

    Copied from Kiss FFT, converted from C++ and modified.

    Copyright (c) 2003-2010, Mark Borgerding

    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
    * Neither the author nor the names of any contributors may be used to
      endorse or promote products derived from this software without specific
      prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """
    if type(n) != int:
        raise ValueError('n must be an integer.')
    if n <= 0:
        raise ValueError('n must be greater than 0.')
    if n & 1 == 1:
        n += 1
    while True:
        m = n
        while m & 1 == 0: m /= 2
        while m % 3 == 0: m /= 3
        while m % 5 == 0: m /= 5
        if m == 1:
            break # n is completely factorable by twos, threes, and fives and is even
        n += 2
    return n
