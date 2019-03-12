#!/bin/sh
# This file is in the public domain.

set -e

for d in ../project/*; do
  mkdir -pv ${d}/output
done
