#!/usr/bin/env Rscript 
#
# Add sample size from .details file to summary data - moba qc
# IB aug 25 

rm(list = ls())

# Command-line parsing 
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript 04_get_sample_size.R <summary_folder>")
}
summary_folder <- args[1]
#details_file   <- args[2]

# Packages 
check_and_install <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
invisible(lapply(c("dplyr","data.table","tools","readr"), check_and_install))

# Helper to read in summary data 
read_data_helper <- function(file) {
  dat <- read_delim(file, col_names = TRUE, na = c("", "NA"),
                    trim_ws = TRUE, show_col_types = FALSE)
  probs <- readr::problems(dat)
  if (nrow(probs) > 0) warning(sprintf("Parsing issues for file %s!", file))
  if (is.null(dat) || nrow(dat) == 0) stop("File is empty or not readable.")
  # convert all except SNP & allele cols to numeric
  num_cols <- setdiff(names(dat), c("Predictor", "A1", "A2"))
  dat[num_cols] <- lapply(dat[num_cols], function(x) suppressWarnings(as.numeric(x)))
  setDT(dat)
  dat
}

# Add sample size column 
add_sample_size_to_summary <- function(summary_folder) {
  
  # List summary files in the summary_folder
  summary_files <- list.files(summary_folder, full.names = TRUE)
  summary_files <- summary_files[!grepl("\\.details$", summary_files)]
  
  if (length(summary_files) == 0) {
    stop("No summary data files found in ", summary_folder)
  }
  
  # Extract Num_Samples from .details file
  details_file <- list.files(summary_folder, pattern = "\\.details$", full.names = TRUE)
  details_lines <- readLines(details_file)
  sample_line   <- details_lines[grepl("^Num_Samples", details_lines, ignore.case = TRUE)]
  if (length(sample_line) == 0) stop("No 'Num_Samples' line found in ", details_file)
  num_samples   <- as.numeric(strsplit(sample_line, "\\s+")[[1]][2])
  if (is.na(num_samples)) stop("Could not parse Num_Samples in ", details_file)
  cat("Extracted Num_Samples:", num_samples, "\n")
  
  # Loop over summary files and add Num_Samples
  for (sum_file in summary_files) {
    cat("Processing:", sum_file, "\n")
    sum_data <- read_data_helper(sum_file)
    sum_data$Num_Samples <- num_samples
    fwrite(sum_data, file = sum_file, sep = "\t", quote = FALSE, na = "NA")
    cat("Updated summary file saved:", sum_file, "\n")
  }
  
  cat("Sample size addition complete.\n")
}

# ---------- 3. Apply function --------------------
add_sample_size_to_summary(summary_folder)