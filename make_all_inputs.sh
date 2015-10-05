#!/usr/bin/env bash

# make_all_inputs.sh: xxx

## This script should be copied into the directory where all of the
## input files are going to be generated!

XYZ_PATH=NewXYZFiles

# What kind of $molecule section do we need?
#  joptype = freq          => supersystem
#  jobtype = freq (w/o CT) => fragments
#  jobtype = eda           => fragments
#  jobtype = eda (w/ COVP) => fragments_covp

declare -A jobtypes_molecule_sections

jobtypes_molecule_sections[freq]=supersystem
jobtypes_molecule_sections[freq_noCT]=fragments
jobtypes_molecule_sections[eda]=fragments
jobtypes_molecule_sections[eda_covp]=fragments_covp

declare -A n_qm_arr

n_qm_arr[supersystem]=$(seq 0 3)
n_qm_arr[fragments]=$(seq 1 3)
n_qm_arr[fragments_covp]=$(seq 1 3)

declare -A python_flags

python_flags[supersystem]="--make-supersystem"
python_flags[fragments]=""
python_flags[fragments_covp]="--qchem-covp"

n_mm_1=$(seq 0 2 16)
n_mm_2=(32 64 128)
n_mm_arr=( ${n_mm_1[@]} ${n_mm_2[@]} )

# The *root* directory for the MBE code.
MBE_ROOT_DIR="${HOME}/development/mbe"
DROPLET_DIR="${MBE_ROOT_DIR}/examples/droplet"

# These files contain the point charges that will be used for the
# cation and anion; the format is such that they can be copied
# directly from a Q-Chem output (Mulliken, Hirshfeld, ChElPG,
# Merz-Kollman...).
pc_output_anion="point_charges_anion.txt"
pc_output_cation="point_charges_cation.txt"

# Make inputs for droplet (CO2 in ionic
# liquid) calculations, varying the number of ionic liquid pairs that
# are treated using QM and MM.
make_input_sections() {
    MOLECULE_SECTION_TYPE="${1}"

    INPUT_SECTIONS_DIR="qchem_molecule_external_charges_stubs_${MOLECULE_SECTION_TYPE}"

    mkdir "${INPUT_SECTIONS_DIR}" || return
    echo "${MOLECULE_SECTION_TYPE}"
    cd "${INPUT_SECTIONS_DIR}"
    cp "${DROPLET_DIR}/${pc_output_anion}" .
    cp "${DROPLET_DIR}/${pc_output_cation}" .
    for n_qm in ${n_qm_arr[${MOLECULE_SECTION_TYPE}]}; do
        for n_mm in ${n_mm_arr[@]}; do
            # echo "${n_qm}" "${n_mm}"
            python_string="python ${DROPLET_DIR}/droplet.py --write-input-sections-qchem --num-closest-pairs-qm=${n_qm} --num-closest-pairs-mm=${n_mm} --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --path=../${XYZ_PATH} ${python_flags[${MOLECULE_SECTION_TYPE}]}"
            echo "${python_string}"
            eval "${python_string}"
        done
        # echo "${n_qm}" "all"
        python_string="python ${DROPLET_DIR}/droplet.py --write-input-sections-qchem --num-closest-pairs-qm=${n_qm} --all-other-pairs-mm --point-charge-output-cation=${pc_output_cation} --point-charge-output-anion=${pc_output_anion} --path=../${XYZ_PATH} ${python_flags[${MOLECULE_SECTION_TYPE}]}"
        echo "${python_string}"
        eval "${python_string}"
        section_dir="sections_${n_qm}qm"
        mkdir "${section_dir}"
        mv ./drop_*_"${n_qm}"qm_* "${section_dir}"
    done
    rm "./${pc_output_anion}"
    rm "./${pc_output_cation}"
    cd ../
}

# Create input files by appending the $rem section to files that
# already contain $molecule/$external_charges sections.
make_inputs() {
    CALC_TYPE="${1}"
    MOLECULE_SECTION_TYPE="${jobtypes_molecule_sections[${CALC_TYPE}]}"

    INPUTS_DIR="inputs_${CALC_TYPE}"
    remfile="droplet_qchem_rem_section_${CALC_TYPE}"

    # So our work isn't overwritten on accident.
    mkdir "${INPUTS_DIR}" || return
    echo "${CALC_TYPE} ${MOLECULE_SECTION_TYPE}"
    cd "${INPUTS_DIR}"

    cp "${DROPLET_DIR}/${remfile}" .
    for n_qm in ${n_qm_arr[${MOLECULE_SECTION_TYPE}]}; do
        n_qm_dir="${INPUTS_DIR}_${n_qm}qm"
        mkdir "${n_qm_dir}"
        cd "${n_qm_dir}"
        echo "${PWD}"
        genfiles="$(find ${PWD}/../../qchem_molecule_external_charges_stubs_${MOLECULE_SECTION_TYPE} -type f -name "drop_*_${n_qm}qm_*mm")"
        for genfile in ${genfiles[@]}; do
            genfilebase=$(basename "${genfile}")
            inputfile="${genfilebase//mm/mm_${CALC_TYPE}.in}"
            cp "${genfile}" "${inputfile}"
            cat "../${remfile}" >> "${inputfile}"
        done
        cd ../
    done
    cd ../
}

##########

cd "${XYZ_PATH}"
eval "python ${DROPLET_DIR}/droplet.py --rename-new"
cd ..

for molecule_section_type in ${jobtypes_molecule_sections[@]}; do
    make_input_sections "${molecule_section_type}"
done

for calc_type in ${!jobtypes_molecule_sections[@]}; do
    make_inputs "${calc_type}"
done
