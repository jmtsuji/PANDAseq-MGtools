#!/bin/bash
# Created April 21, 2016, by Jackson Tsuji
# Description: Uses fasta_length_filter_v0.0.1.R to pull desired lengths from
# Last updated: April 29 2016

# Basic script stuff (from Buffalo's Bioinformatics Data Skills book):
set -e
set -u
set -o pipefail

script_version=0.2.0
# Completely revamping the script from v0.0.2
# Changed to improve speed compared to v0.1.0

# Input variables:
work_dir=$2
pull_length=$3
pull_number=$4
fasta_fileName=$1

echo "Running $(basename $0), version $script_version. **Assumes input FastA file is gzipped."

# Set filename variable
fasta_fileName_base=($(basename $fasta_fileName .fasta.gz))

echo "Pulling reads of desired length >= ${pull_length} from file {fasta_fileName} and subsetting ${pull_number} reads..."
echo ""

# Gunzip FastA file
cd $work_dir
gunzip $fasta_fileName

# Gets all sequences greater than or equal to specified length
# How it works: for all sequence lines of length greater than or equal to pull_length, prints match (i.e. sequence line) and preceeding line before match (i.e. identifier).
# Template awk code from http://www.unixcl.com/2008/05/print-currentnextprevious-line-using.html, accessed April 22, 2016. Idea for -v option to input variable from http://stackoverflow.com/questions/19075671/how-to-use-shell-variables-in-awk-script (by Jotne, edited by Chad), accessed April 22, 2016. AND logical conditional from http://unix.stackexchange.com/questions/216142/how-can-i-use-multiple-if-statement-inside-another-if-statement-of-a-awk-program (by Peter.O), accessed April 22, 2016
awk -v num=${pull_length} ' { if ( !/^>/ && length($0)>=num ) {print x; print};{x=$0} } ' "${fasta_fileName_base}.fasta" > "${fasta_fileName_base}_lengthGEq${pull_length}.fasta"

# Gzip original file
gzip ${fasta_fileName_base}.fasta

# Get number of sequences pulled and print (works with echo text before the start of the loop)
seqNum=$(grep -c "^>" "${fasta_fileName_base}_lengthGEq${pull_length}.fasta")
echo "Sequences pulled of length >= ${pull_length} in ${fasta_fileName_base}_lengthGEq${pull_length}.fasta: ${seqNum}"

# Grab a subset of [pull_number] sequences for indel analysis
seqtk sample -s 57 "${fasta_fileName_base}_lengthGEq${pull_length}.fasta" $pull_number > "${fasta_fileName_base}_lengthGEq${pull_length}_subset${pull_number}.fasta"

echo ""
echo "$(basename $0): Finished."
