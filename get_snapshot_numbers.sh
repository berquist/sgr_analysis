#!/usr/bin/env sh

find . -type f -print0 | xargs -0 -n 1 basename | cut -d "_" -f 2 | sort | uniq
