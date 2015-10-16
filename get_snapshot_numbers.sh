#!/usr/bin/env sh

# get_snapshot_numbers.sh: run inside of a directory with calculation
# inputs to get a sorted list of snapshot numbers

find . -type f -name "*drop*.in" -print0 | \
    xargs -0 -n 1 basename | \
    cut -d "_" -f 2 | \
    sort | \
    uniq
