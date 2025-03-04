#!/bin/bash

# Regrid each component of the basis function decomposition from 1-degree to 2x2.5 degrees

set -aex

source /etc/profile.d/conda.sh
conda activate ./.conda_env

grid_file=$1
IFS=' ' read -ra input_files <<< "$2"
output_directory=$3

remap () {
    local input_file=$1
    local basename=$(basename "$(echo $input_file | grep -oP '.*(?=-\d+\.nc)')")
    local year=$(echo $input_file | grep -oP '(?<=-)\d+(?=.nc)')
    local output_file="${output_directory}/${basename}-2x25-${year}.nc"
    echo "Remapping $input_file to $output_file"
    cdo -w -f nc2 -z zip_6 -remapcon,$grid_file $input_file $output_file
}

# Remap each file in parallel
for file in "${input_files[@]}"; do
    remap $file &
done
wait