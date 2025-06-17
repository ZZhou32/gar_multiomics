#!/bin/bash

# ========================================
# ATAC-seq Peak Counts Matrix Generation
# ========================================

# This script provides multiple methods to generate peak count matrices
# Choose the method that matches your starting data

echo "=== ATAC-seq Peak Counts Matrix Generation ==="

# ========================================
# METHOD 1: From BAM files using bedtools
# ========================================

echo "Method 1: Generating counts from BAM files and peak BED file"

# Prerequisites:
# - Aligned BAM files for each sample
# - Peak regions in BED format (from MACS2 or similar)
# - bedtools installed

# Step 1: Create a merged peak set
echo "Creating merged peak set..."

# If you have individual peak files for each sample:
# cat sample1_peaks.bed sample2_peaks.bed ... > all_peaks.bed
# sort -k1,1 -k2,2n all_peaks.bed > all_peaks_sorted.bed
# bedtools merge -i all_peaks_sorted.bed > merged_peaks.bed

# Or if you have a consensus peak set already:
# cp consensus_peaks.bed merged_peaks.bed

# Step 2: Count reads in peaks for each sample
echo "Counting reads in peaks for each sample..."

PEAK_FILE="merged_peaks.bed"
OUTPUT_DIR="peak_counts"
mkdir -p ${OUTPUT_DIR}

# List of your BAM files (update with your actual file paths)
BAM_FILES=(
  "st_21_22_1_.bam"
  "st_21_22_2_.bam" 
  "st_21_22_3_.bam"
  "st_22_23_1_.bam"
  "st_22_23_2_.bam"
  "st_22_23_3_.bam"
  "st_24_25_1_.bam"
  "st_24_25_2_.bam"
  "st_24_25_3_.bam"
  "st_25_25_1_.bam"
  "st_25_25_2_.bam"
  "st_25_25_3_.bam"
  "st_27_28_1_.bam"
  "st_27_28_2_.bam"
  "st_27_28_3_.bam"
  "st_29_29_1_.bam"
  "st_29_29_2_.bam"
  "st_29_29_3_.bam"
  "st_30_31_1_.bam"
  "st_30_31_2_.bam"
  "st_30_31_3_.bam"
  "st_31_32_1_.bam"
  "st_31_32_2_.bam"
  "st_31_32_3_.bam"
  "st_32_33_1_.bam"
  "st_32_33_2_.bam"
  "st_32_33_3_.bam"
  "st_33_34_1_.bam"
  "st_33_34_2_.bam"
  "st_33_34_3_.bam"
)

# Count reads for each sample
for BAM in "${BAM_FILES[@]}"; do
SAMPLE_NAME=$(basename ${BAM} .bam)
echo "Processing ${SAMPLE_NAME}..."

# Count reads overlapping peaks
bedtools coverage -a ${PEAK_FILE} -b ${BAM} -counts > ${OUTPUT_DIR}/${SAMPLE_NAME}_counts.bed

# Extract just the count column
cut -f4 ${OUTPUT_DIR}/${SAMPLE_NAME}_counts.bed > ${OUTPUT_DIR}/${SAMPLE_NAME}_counts.txt
done

# Step 3: Create peak names
echo "Creating peak identifiers..."
awk '{print $1":"$2"-"$3}' ${PEAK_FILE} > ${OUTPUT_DIR}/peak_names.txt

# Step 4: Combine into matrix
echo "Creating count matrix..."
paste ${OUTPUT_DIR}/peak_names.txt ${OUTPUT_DIR}/*_counts.txt > ${OUTPUT_DIR}/peak_counts_matrix_temp.txt

# Add header
SAMPLE_NAMES=$(for BAM in "${BAM_FILES[@]}"; do basename ${BAM} .bam; done | tr '\n' '\t')
echo -e "peak_id\t${SAMPLE_NAMES}" > ${OUTPUT_DIR}/peak_counts_matrix.txt
cat ${OUTPUT_DIR}/peak_counts_matrix_temp.txt >> ${OUTPUT_DIR}/peak_counts_matrix.txt

echo "Peak counts matrix created: ${OUTPUT_DIR}/peak_counts_matrix.txt"

# ========================================
# METHOD 2: Using featureCounts (more robust)
# ========================================

echo -e "\nMethod 2: Using featureCounts"

# Convert BED to SAF format for featureCounts
echo "Converting BED to SAF format..."
awk 'BEGIN{OFS="\t"; print "GeneID\tChr\tStart\tEnd\tStrand"} 
     {print $1":"$2"-"$3, $1, $2+1, $3, "."}' ${PEAK_FILE} > peaks.saf

# Run featureCounts
echo "Running featureCounts..."
featureCounts \
-F SAF \
-a peaks.saf \
-o peak_counts_featureCounts.txt \
-T 8 \
--minOverlap 1 \
--fracOverlap 0 \
-p \
"${BAM_FILES[@]}"

# Clean up the output
echo "Cleaning featureCounts output..."
# Remove comment lines and extract count columns
sed '1d' peak_counts_featureCounts.txt | cut -f1,7- > gar_atac_peak_counts.txt

echo "Peak counts matrix created: gar_atac_peak_counts.txt"

# ========================================
# METHOD 3: Using DiffBind (R-based)
# ========================================

echo -e "\nMethod 3: R script using DiffBind"

cat > generate_peak_counts_diffbind.R << 'EOF'
#!/usr/bin/env Rscript

# DiffBind method for peak counting
library(DiffBind)
library(GenomicRanges)

# Create sample sheet
samples <- data.frame(
  SampleID = c("st_21_22_1_", "st_21_22_2_", "st_21_22_3_",
               "st_22_23_1_", "st_22_23_2_", "st_22_23_3_",
               "st_24_25_1_", "st_24_25_2_", "st_24_25_3_",
               "st_25_25_1_", "st_25_25_2_", "st_25_25_3_",
               "st_27_28_1_", "st_27_28_2_", "st_27_28_3_",
               "st_29_29_1_", "st_29_29_2_", "st_29_29_3_",
               "st_30_31_1_", "st_30_31_2_", "st_30_31_3_",
               "st_31_32_1_", "st_31_32_2_", "st_31_32_3_",
               "st_32_33_1_", "st_32_33_2_", "st_32_33_3_",
               "st_33_34_1_", "st_33_34_2_", "st_33_34_3_"),
  Tissue = "GAR_tissue",
  Factor = "ATAC",
  Condition = rep(c("21_22", "22_23", "24_25", "25_25", "27_28", 
                    "29_29", "30_31", "31_32", "32_33", "33_34"), each = 3),
  Treatment = rep(c("21_22", "22_23", "24_25", "25_25", "27_28", 
                    "29_29", "30_31", "31_32", "32_33", "33_34"), each = 3),
  Replicate = rep(1:3, 10),
  bamReads = paste0(c("st_21_22_1_", "st_21_22_2_", "st_21_22_3_",
                      "st_22_23_1_", "st_22_23_2_", "st_22_23_3_",
                      "st_24_25_1_", "st_24_25_2_", "st_24_25_3_",
                      "st_25_25_1_", "st_25_25_2_", "st_25_25_3_",
                      "st_27_28_1_", "st_27_28_2_", "st_27_28_3_",
                      "st_29_29_1_", "st_29_29_2_", "st_29_29_3_",
                      "st_30_31_1_", "st_30_31_2_", "st_30_31_3_",
                      "st_31_32_1_", "st_31_32_2_", "st_31_32_3_",
                      "st_32_33_1_", "st_32_33_2_", "st_32_33_3_",
                      "st_33_34_1_", "st_33_34_2_", "st_33_34_3_"), ".bam"),
  Peaks = "consensus_peaks.bed",  # Single consensus peak file
  PeakCaller = "bed"
)

# Write sample sheet
write.csv(samples, "gar_samples.csv", row.names = FALSE)

# Create DiffBind object
dba_obj <- dba(sampleSheet = "gar_samples.csv")

# Count reads in peaks
dba_obj <- dba.count(dba_obj, summits = FALSE)

# Extract count matrix
count_matrix <- dba.peakset(dba_obj, bRetrieve = TRUE, DataType = DBA_DATA_FRAME)

# Format for output
peak_names <- paste(count_matrix$Chr, paste(count_matrix$Start, count_matrix$End, sep = "-"), sep = ":")
count_data <- count_matrix[, 4:ncol(count_matrix)]
rownames(count_data) <- peak_names

# Save count matrix
write.table(count_data, "gar_atac_peak_counts_diffbind.txt", 
            sep = "\t", quote = FALSE, col.names = NA)

cat("DiffBind peak count matrix saved: gar_atac_peak_counts_diffbind.txt\n")
cat("Dimensions:", nrow(count_data), "peaks x", ncol(count_data), "samples\n")

EOF

echo "R script created: generate_peak_counts_diffbind.R"
echo "Run with: Rscript generate_peak_counts_diffbind.R"

# ========================================
# METHOD 4: Using HOMER (if you have tag directories)
# ========================================

echo -e "\nMethod 4: Using HOMER"

cat > generate_peak_counts_homer.sh << 'EOF'
#!/bin/bash

# HOMER method (if you have HOMER tag directories)

PEAK_FILE="merged_peaks.bed"
OUTPUT_DIR="homer_counts"
mkdir -p ${OUTPUT_DIR}

# List of HOMER tag directories
TAG_DIRS=(
  "st_21_22_1_/"
  "st_21_22_2_/"
  "st_21_22_3_/"
  "st_22_23_1_/"
  "st_22_23_2_/"
  "st_22_23_3_/"
  "st_24_25_1_/"
  "st_24_25_2_/"
  "st_24_25_3_/"
  "st_25_25_1_/"
  "st_25_25_2_/"
  "st_25_25_3_/"
  "st_27_28_1_/"
  "st_27_28_2_/"
  "st_27_28_3_/"
  "st_29_29_1_/"
  "st_29_29_2_/"
  "st_29_29_3_/"
  "st_30_31_1_/"
  "st_30_31_2_/"
  "st_30_31_3_/"
  "st_31_32_1_/"
  "st_31_32_2_/"
  "st_31_32_3_/"
  "st_32_33_1_/"
  "st_32_33_2_/"
  "st_32_33_3_/"
  "st_33_34_1_/"
  "st_33_34_2_/"
  "st_33_34_3_/"
)

# Run annotatePeaks.pl to count tags
annotatePeaks.pl ${PEAK_FILE} none -d "${TAG_DIRS[@]}" -raw > ${OUTPUT_DIR}/homer_counts_raw.txt

# Extract just the count columns
cut -f1,8- ${OUTPUT_DIR}/homer_counts_raw.txt > ${OUTPUT_DIR}/homer_counts_clean.txt

echo "HOMER peak counts created: ${OUTPUT_DIR}/homer_counts_clean.txt"

EOF

chmod +x generate_peak_counts_homer.sh
echo "HOMER script created: generate_peak_counts_homer.sh"

# ========================================
# METHOD 5: Quick R script for existing peak/count files
# ========================================

echo -e "\nMethod 5: R script to combine existing count files"

cat > combine_peak_counts.R << 'EOF'
#!/usr/bin/env Rscript

# If you have individual count files for each sample
library(dplyr)

# List of count files (update paths as needed)
count_files <- c(
  "st_21_22_1_counts.txt", "st_21_22_2_counts.txt", "st_21_22_3_counts.txt",
  "st_22_23_1_counts.txt", "st_22_23_2_counts.txt", "st_22_23_3_counts.txt",
  "st_24_25_1_counts.txt", "st_24_25_2_counts.txt", "st_24_25_3_counts.txt",
  "st_25_25_1_counts.txt", "st_25_25_2_counts.txt", "st_25_25_3_counts.txt",
  "st_27_28_1_counts.txt", "st_27_28_2_counts.txt", "st_27_28_3_counts.txt",
  "st_29_29_1_counts.txt", "st_29_29_2_counts.txt", "st_29_29_3_counts.txt",
  "st_30_31_1_counts.txt", "st_30_31_2_counts.txt", "st_30_31_3_counts.txt",
  "st_31_32_1_counts.txt", "st_31_32_2_counts.txt", "st_31_32_3_counts.txt",
  "st_32_33_1_counts.txt", "st_32_33_2_counts.txt", "st_32_33_3_counts.txt",
  "st_33_34_1_counts.txt", "st_33_34_2_counts.txt", "st_33_34_3_counts.txt"
)

# Read peak regions (assuming BED format)
peaks <- read.table("merged_peaks.bed", header = FALSE, 
                    col.names = c("chr", "start", "end"))
peak_names <- paste(peaks$chr, paste(peaks$start, peaks$end, sep = "-"), sep = ":")

# Initialize count matrix
count_matrix <- data.frame(peak_id = peak_names)

# Read and combine count files
for (i in seq_along(count_files)) {
  if (file.exists(count_files[i])) {
    counts <- read.table(count_files[i], header = FALSE)
    sample_name <- gsub("_counts.txt", "", basename(count_files[i]))
    count_matrix[[sample_name]] <- counts[, 1]
  } else {
    cat("Warning: File not found:", count_files[i], "\n")
  }
}

# Save combined matrix
write.table(count_matrix, "gar_atac_peak_counts_combined.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Combined peak count matrix saved: gar_atac_peak_counts_combined.txt\n")
cat("Dimensions:", nrow(count_matrix), "peaks x", ncol(count_matrix)-1, "samples\n")

EOF

echo "R combination script created: combine_peak_counts.R"

# ========================================
# Summary and next steps
# ========================================

echo -e "\n=== Summary ==="
echo "Created scripts for 5 different methods to generate peak count matrices:"
echo "1. bedtools coverage - Basic read counting"
echo "2. featureCounts - Robust read counting (recommended)"
echo "3. DiffBind (R) - Comprehensive ATAC-seq analysis"
echo "4. HOMER - If you have HOMER tag directories"
echo "5. Combine existing files (R) - If you already have individual count files"

echo -e "\nRecommended workflow:"
echo "1. Use Method 2 (featureCounts) for most robust results"
echo "2. Ensure you have:"
echo "   - Aligned BAM files for all 30 samples"
echo "   - Consensus peak set (BED format)"
echo "   - Proper file naming (st_[stage]_[replicate]_.bam)"

echo -e "\nTo run:"
echo "bash generate_peak_counts.sh    # This script"
echo "# OR"
echo "Rscript generate_peak_counts_diffbind.R"

echo -e "\nOutput will be: gar_atac_peak_counts.txt"
echo "Ready for use in the correlation analysis script!"
EOF