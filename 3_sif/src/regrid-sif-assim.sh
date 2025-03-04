#!/bin/bash

# Combine 1-degree hourly SIF and ASSIM files and regrid to 2x2.5 degrees

set -aex

source /etc/profile.d/conda.sh
conda activate ./.conda_env

grid_file=$1
IFS=' ' read -ra assim_input_files <<< "$2"
IFS=' ' read -ra sif_input_files <<< "$3"
output_directory=$4

if [ ${#assim_input_files[@]} -ne ${#sif_input_files[@]} ]; then
    echo "Input lists for SIF and ASSIM files must be the same length."
    exit 1
fi

for i in "${!assim_input_files[@]}"; do

    assim_input_file=${assim_input_files[$i]}
    sif_input_file=${sif_input_files[$i]}

    assim_year=$(echo $assim_input_file | grep -oP '(?<=-)\d+(?=.nc)')
    sif_year=$(echo $sif_input_file | grep -oP '(?<=-)\d+(?=.nc)')

    if [ $sif_year -ne $assim_year ]; then
        echo "SIF and ASSIM yearly files must be aligned by year."
        exit 1
    fi
    echo $assim_year

    output_file="${output_directory}/sib4-hourly-sif-assim-2x25-${assim_year}.nc"

    cdo -w -f nc2 -z zip_6 \
        -remapcon,$grid_file \
        -merge \
        -selvar,sif \
        $sif_input_file \
        -selvar,assim \
        $assim_input_file \
        $output_file

done