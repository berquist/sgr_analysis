#!/usr/bin/env bash

n_qm_arr=$(seq 3)

rm representative_snapshots_*

for n_qm in ${n_qm_arr[@]}; do
    n_qm_dir="inputs_freq_${n_qm}qm_random"
    filename_combined="representative_snapshots_${n_qm}qm"
    echo "${filename_combined}"
    for bindir in "${n_qm_dir}"/*; do
        nbin=$(echo "${bindir}" | cut -d "_" -f 5)
        filename_individual="representative_snapshots_${n_qm}qm_bin_${nbin}"
        snapnums=$(find "${bindir}" -type f -print0 | xargs -0 -n 1 basename | cut -d "_" -f 2 | sort | uniq)
        echo "${filename_individual}"
        echo "${snapnums[@]}" > "${filename_individual}"
        echo "# bin ${nbin}" >> "${filename_combined}"
        echo "${snapnums[@]}" >> "${filename_combined}"
    done
done
