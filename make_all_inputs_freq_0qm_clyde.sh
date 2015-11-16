#!/usr/bin/env bash

make_inputs() {

    # It's probably faster to modify the $rem section and append the
    # modified sections to the $molecule/$external_charges sections,
    # rather than run sed on all the already created inputs.

    inputdir="inputs_freq_0qm_clyde"
    remfile=clyde_rem_section

    # So our work isn't overwritten on accident.
    cd "${inputdir}" || return

    genfiles="$(find /home/eric/Chemistry/calc.sgr/droplets/qchem_molecule_external_charges_stubs_supersystem -type f -name "drop_*_0qm_*mm")"

    for genfile in ${genfiles[@]}; do
        genfilebase=$(basename "${genfile}")
        inputfile="${genfilebase//mm/mm_freq.in}"
        cp "${genfile}" "${inputfile}"
        cat "${remfile}" >> "${inputfile}"
    done

    cd ../
}

####################

make_inputs
