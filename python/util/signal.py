
import numpy as np
from numpy.fft import rfft, irfft

SPIKES_MIN_OFFSET_COEF = 0.35
SPIKES_MIN_ABS_OFFSET = 0.03

def filter_spikes(s):
    flag = np.zeros(len(s))
    for i in range(1, len(s) - 1):
        min_offset = abs(s[i]) * SPIKES_MIN_OFFSET_COEF + SPIKES_MIN_ABS_OFFSET
        if s[i] - s[i-1] > min_offset and s[i] - s[i+1] > min_offset:
            flag[i] = 1.0
        elif s[i-1] - s[i] > min_offset and s[i+1] - s[i] > min_offset:
            flag[i] = -1.0
    s_out = s.copy()
    for i in range(1, len(s_out) - 1):
        if flag[i] != 0.0:
            if flag[i-1] == 0.0 and flag[i+1] == 0.0:
                s_out[i] = 0.5 * (s_out[i-1] + s_out[i+1])
            else:
                s_out[i] = 0.0
    return s_out

def sinc_lowpass(s, fs, cutoff):
    n = len(s)
    if n % 2 == 1:
        n = n - 1
    f = rfft(s, n) # f size: (n/2)+1
    nc = int((cutoff / fs) * n)
    if nc >= (n // 2) + 1:
        raise ValueError("Cutoff frequency is too high.")
    f[nc:] = 0.0
    return irfft(f)
