#!/bin/bash
# Created April 7, 2016, by Jackson Tsuji
# Description: prints names of sequences in a FastQ file
# Last updated: April 17, 2016

# Basic script stuff (from Buffalo's Bioinformatics Data Skills book):
set -e
set -u
set -o pipefail

script_version=0.0.2
# This version is not verbose but prints FastQ names to STDOUT.

# Processing user input variable
input_file=$1
#output_file=$2

#echo "Running $(basename $0), version $script_version, on input FastQ file $input_file."

awk '{ if (NR%4 == 1) { print $0 } }' $input_file | cut -f 2- -d '@'
# awk code is using the "if, then" format: IF ? THEN (: ELSE)
# General idea: if row number is a multiple of four, then rename to the "id" variable appended by a numerical counter. Otherwise (i.e. row is not a multiple of four), print the entire original line


# Some references:
##  Main awk code adapted from https://www.biostars.org/p/68477/ (accessed April 7, 2016), from Frederic Mahe

#echo "$(basename $0): Finished."
