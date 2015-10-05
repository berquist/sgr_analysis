#!/usr/bin/env bash

methods=(
    b3lyp
    blyp
    tpss
    hf
    wb97x-d
    ri-mp2
)

declare -A basis_sets

basis_sets[6-31gdp]="6-31g**"
basis_sets[cc-pvtz]="cc-pvtz"

# Newlines should separate each entry, but a newline is unnecessary
# after the last entry.
new_rem_lines=" ideriv = 1\n aux_basis = rimp2-cc-pvqz\n n_frozen_core = fc"

make_method_basis_set_inputs() {
    method="${1}"
    basis_set="${2}"

    # It's probably faster to modify the $rem section and append the
    # modified sections to the $molecule/$external_charges sections,
    # rather than run sed on all the already created inputs.

    inputdir="inputs_freq_0qm_${method}_${basis_set}"
    orig_remfile=droplet_qchem_rem_section_freq
    remfile="${orig_remfile}_${method}_${basis_set}"

    # So our work isn't overwritten on accident.
    mkdir "${inputdir}" || return
    echo "${inputdir}"
    cd "${inputdir}"

    cp ../"${orig_remfile}" "${remfile}"
    sed -i \
        -e "s/\$rem/\$rem\n${new_rem_lines}/" \
        -e "s/method = b3lyp/method = ${method}/" \
        -e "s/basis = 6-31g\*\*/basis = ${basis_sets[${basis_set}]}/" \
        "${remfile}"

    genfiles="$(find "${PWD}"/../../qchem_molecule_external_charges_stubs_supersystem -type f -name "drop_*_0qm_*mm")"

    for genfile in ${genfiles[@]}; do
        genfilebase=$(basename "${genfile}")
        inputfile="${genfilebase//mm/mm_freq.in}"
        cp "${genfile}" "${inputfile}"
        cat "${remfile}" >> "${inputfile}"
    done

    cd ../
}

####################

cd inputs_freq

for method in ${methods[@]}; do
    for basis_set in ${!basis_sets[@]}; do
        make_method_basis_set_inputs "${method}" "${basis_set}"
    done
done
