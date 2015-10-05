#!/usr/bin/env bash

for f in $(ls .); do
    echo "${f}"
    tar cf "${f}.tar" "${f}"
    echo "${f}.tar"
    pigz --best -p 4 "${f}.tar"
    rm -r "${f}"
done
