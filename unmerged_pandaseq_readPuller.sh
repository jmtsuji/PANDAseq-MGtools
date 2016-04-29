#!/bin/bash
# Created April 17, 2016, by Jackson Tsuji
# Description: uses a collection of scripts to quickly pull unaligned PANDAseq reads from original metagenome read files.
# Last updated: April 29, 2016

# Basic script stuff (from Buffalo's Bioinformatics Data Skills book):
set -e
set -u
set -o pipefail

script_version=0.0.1

# Input variables:
work_folder=$4
aligned_pandaseq_FastA_fileName=$1
R1_rawReadFile_name=$2
R2_rawReadFile_name=$3

echo "Running $(basename $0), version $script_version. **NOTE: assumes PANDAseq aligned reads are gzipped in FastA format and that raw reads are gzipped in FastQ format."

## Preamble: get base names of input files for use within the script
aligned_pandaseq_FastA_fileName_base=($(basename $aligned_pandaseq_FastA_fileName .fasta.gz))
R1_rawReadFile_name_base=($(basename $R1_rawReadFile_name .fastq.gz))
R2_rawReadFile_name_base=($(basename $R2_rawReadFile_name .fastq.gz))

## Preamble: unzip source files
echo "Unzipping source files into working directory..."
mkdir -p $work_folder
cd $work_folder
gunzip "${work_folder}/${aligned_pandaseq_FastA_fileName}"
gunzip "${work_folder}/${R1_rawReadFile_name}"
gunzip "${work_folder}/${R2_rawReadFile_name}"

## Step 1: create a list of unaligned sequences to extract from the original raw read files
# Get names of all sequence identifiers for raw reads (using R1 file, arbitrarily), renamed to standard identifier
echo "Determining sequence identifiers of unaligned reads to prepare for extraction..."
fastq_getNames_v0.0.2.sh "${R1_rawReadFile_name_base}.fastq" | cut -f 1 -d " " > "${R1_rawReadFile_name_base}_names.txt"
# Get names of PANDAseq sequences, renamed to standard identifier
grep "^>" "${aligned_pandaseq_FastA_fileName_base}.fasta" | cut -f 2 -d ">" | cut -f 1-7 -d ":" > "${aligned_pandaseq_FastA_fileName_base}_names.txt"
# Find unaligned sequences (i.e. IDs not in the aligned file names):
awk 'FNR==NR {a[$0]++; next} !a[$0]' "${aligned_pandaseq_FastA_fileName_base}_names.txt" "${R1_rawReadFile_name_base}_names.txt" > "${R1_rawReadFile_name_base}_unmerged_names.txt"
# Citation for awk code: suggested on http://stackoverflow.com/questions/4717250/extracting-unique-values-between-2-sets-files (accessed April 13, 2016) by SiegeX, edited by Rorick

echo ""
echo "IDs of unaligned reads have been found. Check to make sure the numbers add up (first two should sum to the third):"
wc -l "${aligned_pandaseq_FastA_fileName_base}_names.txt"
wc -l "${R1_rawReadFile_name_base}_unmerged_names.txt"
wc -l "${R1_rawReadFile_name_base}_names.txt"

## Step 2: Pull unaligned reads out of raw read files
echo ""
echo "Using unaligned read ID list to pull reads out of R1 and R2 read files, using seqtk subseq command..."
seqtk subseq "${R1_rawReadFile_name_base}.fastq" "${R1_rawReadFile_name_base}_unmerged_names.txt" > "${R1_rawReadFile_name_base}_unmerged.fastq"
seqtk subseq "${R2_rawReadFile_name_base}.fastq" "${R1_rawReadFile_name_base}_unmerged_names.txt" > "${R2_rawReadFile_name_base}_unmerged.fastq"
# Side note: it's okay that the "names" file has the R1 identifier in the R2 extraction code above; names should be the same between the R1 and R2 files for extraction.

echo ""

# Cleanup
echo "Cleanup: compressing source and output files... also compressing sequence ID lists for future reference..."
gzip "${R1_rawReadFile_name_base}_unmerged.fastq" "${R2_rawReadFile_name_base}_unmerged.fastq"
gzip "${aligned_pandaseq_FastA_fileName_base}.fasta" "${R1_rawReadFile_name_base}.fastq" "${R2_rawReadFile_name_base}.fastq"
gzip "${R1_rawReadFile_name_base}_names.txt" "${aligned_pandaseq_FastA_fileName_base}_names.txt" "${R1_rawReadFile_name_base}_unmerged_names.txt"
echo ""
echo ""

echo "$(basename $0): Finished."
