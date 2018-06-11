#!/bin/bash
# This file is in the public domain.

set -e

rsh_cmd="ssh -l ultrassom"
server=bolha

rsync -vaz --progress --delete-before --stats -e "${rsh_cmd}" \
  --exclude='.git*' --include='*/' --include='*.h5' --exclude='*' \
  ${server}:/home/ultrassom/project/qtcreator/us_lab4a/project .

sync

echo === TRANSFER OK
