#!/bin/bash

pesmatdirs=(PESmats RigidPESmats)

for pesmatdir in "${pesmatdirs[@]}"; do
    cd "${pesmatdir}" || exit
    ls *.dat | xargs python ../../analysis/unused/plot_pesmats.py
    cd .. || exit
done
