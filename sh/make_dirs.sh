#!/bin/sh

set -e

for d in ../project/*; do
  mkdir -p ${d}/output
done
