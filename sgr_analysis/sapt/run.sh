#!/bin/bash

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

directories=(
    './631gdp/ct/'
    './jun-cc-pvdz/'
    './jun-cc-pvtz/ct/'
)

SAPT_CSV=data_sapt.csv
SAPT_XLS=data_sapt.xlsx

rm ${SCRIPTDIR}/${SAPT_CSV}
rm ${SCRIPTDIR}/${SAPT_XLS}

for directory in ${directories[@]}; do
    pushd ${directory}
    python ${SCRIPTDIR}/collect_sapt_data.py --averages ./*.out
    cat ${SAPT_CSV} >> ${SCRIPTDIR}/${SAPT_CSV}
    # python concat_xlsx.py --output=${SCRIPTDIR}/${SAPT_XLS} ${SAPT_XLS}
    popd
done
