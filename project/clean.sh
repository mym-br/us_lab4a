#!/bin/bash
# This file is in the public domain.

set -e

find . -type d -name 'output*' -prune -exec echo '{}' \; -exec rm -rI '{}' \;
