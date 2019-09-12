#!/bin/bash
# This file is in the public domain.

set -e

find . -type d -name 'output*' -exec rm -rf '{}' \;
