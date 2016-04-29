# PANDAseq-MGtools
Collection of scripts to optimize and use PANDAseq for paired-end merging of metagenomic reads

Copyright Jackson Tsuji and Katherine Naismith, 2016

*__NOTE__: All scripts are currently in early development and are not guaranteed to function properly. Please contact the author before use.

These scripts are designed to both optimize PANDAseq paired-end read merging for a specific metagenome dataset of interest and to faciliate the use of PANDAseq output for metagenomic assembly.

## Dependencies
- PANDAseq (recommended verson 2.8 or higher)
- seqtk
- R
- Python

## Installation
To use these scripts, simply copy them into a desired directory and add that directory to your $PATH. It is recommended to use these scripts on a Linux system.

## Usage
Only four of the collection of scripts require user input. The others are required to support analysis performed by these four.

### 1. Recursive optimization of minimum length overlap threshold (-o parameter)
```
merge_opt_commands.sh [run_name] [fwd_read_filepath.fastq.gz] [rev_read_filepath.fastq.gz] [/output/data/directory] > logfile.txt
```
Runs PANDAseq on the given input read files 25 times, from minimum output threshold length of 1 to 25, using a quality threshold of 0.9 and otherwise default parameters. Number of iterations, as well as quality threshold, can be changed within the script on lines 25 and 28, respectively.
Output: generates plots in the output folder showing the distrubiton of sequence lengths output in the aligned read file in each run, along with CSV files summarizing the length distribution for easy plotting.

### 2. Checking for suspected insertion-deletion errors in PANDAseq output
A. Firstly, pull out sequence reads of interest from the merged pe-read PANDAseq output file. I recommend the longest reads possible (e.g. 400 and higher if merging 201 bp reads with 1 bp minimum overlap), since these are expected to be the most prone to having indel errors.
```
pull_by_length.sh [merged_fasta_filename.fasta.gz] [/output/data/directory] [pull_length] [subset_number] > logfile.txt
```
Output: gives all sequences in the input file of length >= [pull_length], as well as a subset of that file with [subset_number] random sequences that is convenient for use in BLASTN searches. **I recommend a subset number of 1000.

B. Next, upload the subset FastA file to BLASTN (http://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome) and perform a search (recommended to use the nr database, and default parameters should be fine for most purposes). When finished, scroll to the first sequence with a hit, select "Download" (at the top of the page), and downloaded the "Hit Table(text)". This will be used in the next step.

C. Search for indel errors
**__Manually input__ the name of the downloaded hit table at line 54 of the python script "indel_detector.py" (e.g. by opening it in a text editor), as well as the minimum overlap setting you used for this read merge at line 55. I recommend leaving the "buffer" at line 56 set to 10 -- this account for errors in BLASTN alignments extending slightly beyond where they ought to end based on sequence homology.
Then, __move the python script to the folder containing the hit table file__ and run:
```
python indel_detector.py
```
Output: a tab-delimited file summarizing whether each query sequence likely contains an indel ("Gap" or "Overlap") is successfully aligned ("Alignment Succeeded"), or had uninformative search results.

### 3. Retrieving paired-end reads unmerged by PANDAseq for metageomic assembly (note the that the built-in feature in PANDAseq to output these unmerged reads may not function reliably)
After determining optimum settings for PANDAseq, run the following script to retrieve unmerged paired end reads for downstream use:
```
unmerged_pandaseq_readPuller.sh [pandaseq_merged_read_filepath.fasta.gz] [fwd_read_filepath.fastq.gz] [rev_read_filepath.fastq.gz] [/output/data/directory] > logfile.txt
```
