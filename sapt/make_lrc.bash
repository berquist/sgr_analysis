#!/bin/bash

rm -r lrc
mkdir lrc
cp ./*b3lyp_cp.in lrc
cd lrc || exit

declare -A method_dispersion

method_dispersion["b3lyp-d2"]="empirical_grimme"
method_dispersion["b3lyp-d3"]="empirical_grimme3"

declare -A method_basis

method_basis["wb97x-d"]="cc-pvtz"
method_basis["wb97x-d3"]="aug-cc-pvtz"

for inputfile in ./*b3lyp_cp.in; do
    for method in ${!method_dispersion[@]}; do
        newinputfile=${inputfile//b3lyp/${method}_6-31gss}
        cp ${inputfile} ${newinputfile}
        new_rem_lines=" dft_d = ${method_dispersion[${method}]}"
        sed -i "s/\$rem/\$rem\n${new_rem_lines}/" ${newinputfile}
        newinputfile_nocp=${newinputfile//_cp/_nocp}
        cp ${newinputfile} ${newinputfile_nocp}
        sed -i "/eda_bsse/d" ${newinputfile_nocp}
    done
    for method in ${!method_basis[@]}; do
        basis=${method_basis[${method}]}
        newinputfile=${inputfile//b3lyp/${method}_${basis}}
        cp ${inputfile} ${newinputfile}
        sed -i "s/method = b3lyp/method = ${method}/g" ${newinputfile}
        sed -i "s/basis = 6-31g\*\*/basis = ${basis}/g" ${newinputfile}
        newinputfile_nocp=${newinputfile//_cp/_nocp}
        cp ${newinputfile} ${newinputfile_nocp}
        sed -i "/eda_bsse/d" ${newinputfile_nocp}
    done
    rm ${inputfile}
done
