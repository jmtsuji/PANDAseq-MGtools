#!/usr/bin/Rscript --vanilla
# FastA length analyzer
# Description: input a FastA file, and this script outputs the lengths of the sequences, and, if desired, summary stats (NOT yet implemented) and a plot of their distribution.
# Created by Jackson Tsuji on March 18th, 2016
# Note: idea for shebang line is presented at http://stackoverflow.com/questions/12540138/how-can-i-make-my-r-session-vanilla (accessed April 5th, 2016)

inputArgs <- commandArgs(trailingOnly = TRUE)
if (length(inputArgs) < 2) {
  stop("FastA Length Analyzer, version 1.0.9, by Jackson Tsuji (Neufeld Research Group; jackson.tsuji@uwaterloo.ca) \n
        Description: outputs a tabular and graphical summary of the lengths of sequences in a FastA file. Requires R and the ggplot2 R library.
        This script needs at least two input arguments. Must provide (in this order): \n
       1. The input filename (should be a FastA file with multiple sequences; DNA sequence data for each entry much be on a single line rather than split across multiple lines). \n
       2. A unique identifier to append to all output files for this run. \n
       Optional: 3. Output directory path for generated files. *Note that if you specify this, you need to specify the full filepath to the input file for argument 1. \n
       Optional: 4. The maximum sequence length to display in the length distribution plot (although the length table will still contain all entries). \n
       Optional: 5. The minimum sequence length to display in the length distribution plot (although the length table will still contain all entries). ", call. = FALSE)
}
# From http://www.r-bloggers.com/passing-arguments-to-an-r-script-from-command-lines/, accessed approx. March 21, 2016

#############################
# Input variables (can tweak settings here if needed):
fasta_filename <- inputArgs[1]
printRawSeqLengths <- FALSE
output_length_filename <- paste("lengths_", inputArgs[2], ".txt", sep = "")
printLengthTable <- TRUE
output_length_table_filename <- paste("lengths_table_", inputArgs[2], ".csv", sep = "")
printPDF <- TRUE
PDF_name <- paste("lengths_plot_", inputArgs[2], ".pdf", sep = "")
output_dir <- inputArgs[3]
max_length_to_print <- as.numeric(inputArgs[4])
min_length_to_print <- as.numeric(inputArgs[5])
#############################

# Set working directory (where files will be output to) only if supplied by user:
if (is.na(output_dir) == FALSE) {
  setwd(output_dir)
}

# Use bash commands to read in the lengths of all sequences in the FastA file
lengths <- system(paste("grep -v \'^>\' ", fasta_filename, " | awk \'{ print length($0) }\'", sep = ""), intern = TRUE)
# Convert length vector to numeric format for further processing
lengths <- as.numeric(lengths)

# If desired, print a raw list of the lengths of each sequence in the file:
if (printRawSeqLengths == TRUE) {
  write(x = lengths, file = output_length_filename, ncolumns = 1)
}
# summary(lengths)

# If desired, print a tabular summary of the number of sequences of each sequence length
if (printLengthTable == TRUE) {
  # Create length table
  table_lengths <- as.data.frame(table(lengths))
  # Write table to a file
  write.table(x = table_lengths, file = output_length_table_filename, sep = ",", row.names = FALSE, col.names = c("Sequence_length_nt", "Number_of_sequences"))
}

# If desired, print a PDF plot showing the distribution of sequence lengths (*can be a bit slow)
if (printPDF == TRUE) {
  # Preparing for plotting with ggplot2 by converting the raw sequence lengths vector to a dataframe
  lengths_dataframe <- as.data.frame(lengths)
  
  # Specify ggplot title (will be modified if max/min length settings are turned on, in the lines below)
  plot_title <- paste("Sequence length distribution for file:\n", fasta_filename, sep = "")
  
  # If specified by the user, remove lengths above a cutoff value (and also adjust plot title)
  if (is.na(max_length_to_print) == FALSE) {
    lengths_dataframe <- subset(lengths_dataframe, lengths <= max_length_to_print)
    plot_title <- paste(plot_title, "\n", "*Note: only showing sequences under a length of ", max_length_to_print, sep = "")
  }
  
  # If specified by the user, remove lengths below a cutoff value (and also adjust plot title)
  if (is.na(min_length_to_print) == FALSE) {
    lengths_dataframe <- subset(lengths_dataframe, lengths >= min_length_to_print)
    # Determine plot title based on whether or not the title has already been appended above
    if (is.na(max_length_to_print) == FALSE) {
      plot_title <- paste(plot_title, " and over a length of ", min_length_to_print, sep = "")
    }
    else {
      plot_title <- paste(plot_title, "\n", "*Note: only showing sequences over a length of ", min_length_to_print, sep = "")
    }
  }
  
  # Install the ggplot2 library, if needed
  if ("ggplot2" %in% rownames(installed.packages()) == FALSE) {
    install.packages("ggplot2", repos = "http://cran.us.r-project.org")
  }
  # Load the ggplot2 library
  library(ggplot2)#, lib.loc = library_path)
  
  # Create the ggplot
  length_plot <- ggplot(lengths_dataframe, aes(lengths)) +
    geom_bar() +
    theme_bw() +
    theme(title = element_text(size = 8), axis.title = element_text(size = 11), panel.grid = element_blank()) +
    xlab("DNA sequence length (nt)") +
    ylab("Number of sequences") +
    ggtitle(plot_title)
  
  # Save the ggplot to a PDF
  ggsave(PDF_name, width = 9, height = 7)
}
