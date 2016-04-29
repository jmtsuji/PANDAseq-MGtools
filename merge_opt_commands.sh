#!/bin/bash
# Created March 31, 2016, by Jackson Tsuji
# Description: mini-pipeline for alignment test commands to run iteratively
# Last updated: April 8, 2016

# Basic script stuff (from Buffalo's Bioinformatics Data Skills book):
set -e
set -u
set -o pipefail

script_version=0.0.2a

# Making simple output directory paths:
run_name=$1
input_fwd_filepath=$2
input_rev_filepath=$3
output_data_dir=$4

echo "Running $(basename $0), version $script_version."
echo ""

# Check that the correct output directory exists
mkdir -p "${output_data_dir}"

for i in {1..25}
do
    # Set variables
    qual_threshold=0.9
    min_overlap=$i

#    echo "Aligning using PANDAseq with t=${qual_threshold} and o=${min_overlap}..."

    # Run alignTest script
    merge_test_support.sh $run_name $input_fwd_filepath $input_rev_filepath $output_data_dir $qual_threshold $min_overlap
done

echo ""
echo "$(basename $0): finished."
