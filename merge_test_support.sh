#!/bin/bash
# alignment testing script
# Created March 31st, 2016, by Jackson Tsuji
# Description: Part 3b of pipeline for processing L227 and L442 metagenome reads. For use in testing different settings in aligning reads.
# Last updated: April 29, 2016

# Basic script stuff (from Buffalo's Bioinformatics Data Skills book):
set -e
set -u
set -o pipefail

script_version=0.0.2a

# Making input variables more understandable (passed from pipeline.sh):
run_name=$1
input_fwd_filepath=$2
input_rev_filepath=$3
output_data_dir=$4
qual_threshold=$5
min_overlap=$6

#echo "$(basename $0) version $script_version"

#### Step 1: preliminary setup
# Check that the correct output directory exists
mkdir -p "${output_data_dir}"

#echo "Running for sample ${run_name} with t=${qual_threshold} and o=${min_overlap}."

# Aligning with PandaSeq
echo "Aligning using PANDAseq with t=${qual_threshold} and o=${min_overlap}..."
pandaseq -f $input_fwd_filepath -r $input_rev_filepath -A simple_bayesian -G "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_log.txt.bz2" -T 4 -u "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_unpaired.fasta" -w "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_aligned.fasta" -t ${qual_threshold} -l 200 -o ${min_overlap}

# Printing number of reads in output files
#echo ""
echo "PANDAseq output stats with t=${qual_threshold} and o=${min_overlap}:"
num_aligned=$(grep -c "^>" "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_aligned.fasta")
echo "Number of aligned reads: ${num_aligned}"
num_unpaired=$(grep -c "^>" "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_unpaired.fasta")
echo "Number of unpaired reads: ${num_unpaired}"
echo ""

# Generating length distribution report for PANDAseq aligned output
fasta_length_analyze_v1.0.9.R "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_aligned.fasta" "${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_aligned" "${output_data_dir}"

# Compressing Sickle and PandaSeq output files
gzip --fast "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_unpaired.fasta"
gzip --fast "${output_data_dir}/${run_name}_pandaseq_t${qual_threshold}o${min_overlap}_aligned.fasta"

#echo "Finished running $(basename $0) on sample ${run_name} with t=${qual_threshold} and o=${min_overlap}."
