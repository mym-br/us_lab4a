#!/bin/sh

set -e

in_dir=saved_acquisition
out_dir=saved_acquisition-2
dataset=signal
fs=4.0e6
cutoff_freq=800.0e3

rm -rf ${out_dir}

#../../python/filter_spikes.py single ${in_dir} ${out_dir} \
#  0062/signals-base0016.h5 ${dataset} ${fs} ${cutoff_freq}

../../python/filter_spikes.py sequence ${in_dir} ${out_dir} \
  ${dataset} ${fs} ${cutoff_freq}
cp ${in_dir}/time.h5 ${out_dir}
cp ${in_dir}/y.h5    ${out_dir}
