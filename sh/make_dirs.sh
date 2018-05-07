#!/bin/sh

set -e

for d in ../project/*; do
  mkdir -pv ${d}/output
done
