#!/usr/bin/env Rscript

# Load required libraries
library(rtracklayer)
library(readr)
library(dplyr)

# Function to extract gene positions from gar genome
extract_gene_positions <- function(csv_file, gar_gff_file, output_file = "genes_of_interest.gff3") {
  
  cat("Reading CSV file:", csv_file, "\n")
  
  # Read the CSV file
  de_results <- read_csv(csv_file)
  
  # Check if required columns exist
  if (!all(c("status", "gene") %in% colnames(de_results))) {
    stop("CSV file must contain 'status' and 'gene' columns")
  }
  
  # Filter for genes that are not "nonDE"
  genes_of_interest <- de_results %>%
    filter(status != "nonDE") %>%
    pull(gene) %>%
    unique()
  
  cat("Found", length(genes_of_interest), "genes of interest\n")
  cat("Genes:", paste(head(genes_of_interest, 10), collapse = ", "), 
      ifelse(length(genes_of_interest) > 10, "...", ""), "\n")
  
  # Read the gar genome annotation
  cat("Reading gar genome annotation:", gar_gff_file, "\n")
  gar_annotation <- import(gar_gff_file)
  
  # Extract gene positions
  cat("Extracting gene positions...\n")
  
  # Try different possible gene ID fields
  gene_matches <- c()
  
  # Check common gene identifier fields
  if ("gene_id" %in% names(mcols(gar_annotation))) {
    gene_matches <- gar_annotation[gar_annotation$gene_id %in% genes_of_interest]
  }
  
  if (length(gene_matches) == 0 && "Name" %in% names(mcols(gar_annotation))) {
    gene_matches <- gar_annotation[gar_annotation$Name %in% genes_of_interest]
  }
  
  if (length(gene_matches) == 0 && "ID" %in% names(mcols(gar_annotation))) {
    gene_matches <- gar_annotation[gar_annotation$ID %in% genes_of_interest]
  }
  
  # If still no matches, try partial matching on any attribute
  if (length(gene_matches) == 0) {
    cat("Trying partial matching...\n")
    all_attrs <- mcols(gar_annotation)
    for (gene in genes_of_interest) {
      for (col in names(all_attrs)) {
        if (is.character(all_attrs[[col]])) {
          matches <- grep(gene, all_attrs[[col]], fixed = TRUE)
          if (length(matches) > 0) {
            gene_matches <- c(gene_matches, gar_annotation[matches])
          }
        }
      }
    }
    gene_matches <- unique(gene_matches)
  }
  
  cat("Found", length(gene_matches), "matching genes in the annotation\n")
  
  if (length(gene_matches) == 0) {
    cat("No matching genes found. Checking available gene identifiers...\n")
    
    # Show sample of available gene IDs to help troubleshoot
    if ("gene_id" %in% names(mcols(gar_annotation))) {
      sample_ids <- head(unique(gar_annotation$gene_id), 10)
      cat("Sample gene_ids:", paste(sample_ids, collapse = ", "), "\n")
    }
    if ("Name" %in% names(mcols(gar_annotation))) {
      sample_names <- head(unique(gar_annotation$Name), 10)
      cat("Sample Names:", paste(sample_names, collapse = ", "), "\n")
    }
    
    stop("No matching genes found in annotation file")
  }
  
  # Filter for gene features only (optional, remove if you want all features)
  if ("type" %in% names(mcols(gene_matches))) {
    gene_features <- gene_matches[gene_matches$type == "gene"]
    if (length(gene_features) > 0) {
      gene_matches <- gene_features
      cat("Filtered to gene features only:", length(gene_matches), "entries\n")
    }
  }
  
  # Add status information to the annotations
  gene_matches$differential_status <- "DE"
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    cat("Creating output directory:", output_dir, "\n")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Export to GFF3 file
  cat("Exporting to:", output_file, "\n")
  tryCatch({
    export(gene_matches, output_file, format = "gff3")
  }, error = function(e) {
    cat("Error writing to output file. Trying alternative location...\n")
    # Try writing to current working directory as backup
    backup_file <- file.path(getwd(), basename(output_file))
    cat("Writing to backup location:", backup_file, "\n")
    export(gene_matches, backup_file, format = "gff3")
    output_file <<- backup_file  # Update the output_file variable
  })
  
  # Print summary
  cat("\nSummary:\n")
  cat("- Total genes of interest:", length(genes_of_interest), "\n")
  cat("- Genes found in annotation:", length(gene_matches), "\n")
  cat("- Output file:", output_file, "\n")
  
  # Show coordinates of first few genes
  if (length(gene_matches) > 0) {
    cat("\nFirst few gene coordinates:\n")
    for (i in 1:min(5, length(gene_matches))) {
      gene <- gene_matches[i]
      gene_name <- ifelse("gene_id" %in% names(mcols(gene)), gene$gene_id, 
                          ifelse("Name" %in% names(mcols(gene)), gene$Name, 
                                 ifelse("ID" %in% names(mcols(gene)), gene$ID, "unknown")))
      cat(sprintf("%s: %s:%d-%d (%s)\n", 
                  gene_name, 
                  as.character(seqnames(gene)), 
                  start(gene), 
                  end(gene), 
                  as.character(strand(gene))))
    }
  }
  
  return(gene_matches)
}

# ===================================================================
# MODIFY THESE FILE PATHS FOR YOUR ANALYSIS
# ===================================================================

# Input CSV file with DE results (must have 'status' and 'gene' columns)
csv_file <- "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/example_29_30/2929-3031.csv"

# Gar genome annotation file (GFF3 or GTF format)
gar_gff_file <- "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/gar_genome.gff3"

# Output file for genes of interest
output_file <- "/mnt/ufs18/rs-032/FishEvoDevoGeno/raw_reads_4_tony/salmon_quant/genes_of_interest.gff3"

# ===================================================================
# MAIN EXECUTION - NO NEED TO MODIFY BELOW THIS LINE
# ===================================================================

# Check if running from command line or interactively
if (!interactive()) {
  # Command line usage (optional - overrides the paths above)
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) >= 1) {
    csv_file <- args[1]
  }
  if (length(args) >= 2) {
    gar_gff_file <- args[2]
  }
  if (length(args) >= 3) {
    output_file <- args[3]
  }
}

# Display the file paths being used
cat("=== File Paths ===\n")
cat("CSV file:", csv_file, "\n")
cat("GAR annotation:", gar_gff_file, "\n")
cat("Output file:", output_file, "\n")
cat("==================\n\n")

# Check if files exist
if (!file.exists(csv_file)) {
  stop("CSV file not found: ", csv_file)
}

if (!file.exists(gar_gff_file)) {
  cat("GFF file not found:", gar_gff_file, "\n")
  cat("Let's check what files are available in that directory...\n")
  
  # List files in the directory to help find the correct annotation file
  dir_path <- dirname(gar_gff_file)
  if (dir.exists(dir_path)) {
    cat("Files in", dir_path, ":\n")
    files <- list.files(dir_path, pattern = "\\.(gff|gff3|gtf)$", full.names = TRUE)
    if (length(files) > 0) {
      cat("Available annotation files:\n")
      for (f in files) {
        cat(" -", f, "\n")
      }
      cat("\nPlease update the 'gar_gff_file' variable with the correct path.\n")
    } else {
      cat("No .gff, .gff3, or .gtf files found in this directory.\n")
      cat("All files in directory:\n")
      all_files <- list.files(dir_path, full.names = TRUE)
      for (f in head(all_files, 20)) {  # Show first 20 files
        cat(" -", f, "\n")
      }
      if (length(all_files) > 20) {
        cat(" ... and", length(all_files) - 20, "more files\n")
      }
    }
  } else {
    cat("Directory does not exist:", dir_path, "\n")
  }
  stop("Please correct the gar_gff_file path and try again.")
}

# Run the extraction
cat("Starting gene extraction...\n")
result <- extract_gene_positions(csv_file, gar_gff_file, output_file)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:", output_file, "\n")
cat("You can now load this file into IGV!\n")