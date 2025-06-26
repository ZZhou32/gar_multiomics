library(readr)
BiocManager::install("GRaNIE")
library(GRaNIE)

## Getting data from the tutorial.a
# We load the example data directly from the web:
file_peaks = "https://www.embl.de/download/zaugg/GRaNIE/countsATAC.filtered.tsv.gz"
file_RNA = "https://www.embl.de/download/zaugg/GRaNIE/countsRNA.filtered.tsv.gz"
file_sampleMetadata = "https://www.embl.de/download/zaugg/GRaNIE/sampleMetadata.tsv.gz"

countsRNA.df = read_tsv(file_RNA, col_types = cols())
countsPeaks.df = read_tsv(file_peaks, col_types = cols())
sampleMetadata.df = read_tsv(file_sampleMetadata, col_types = cols())


# Let's check how the data looks like
countsRNA.df
countsPeaks.df
sampleMetadata.df

# Save the name of the respective ID columns
idColumn_peaks = "peakID"
idColumn_RNA = "ENSEMBL"

genomeAssembly = "hg38"  #Either hg19, hg38 or mm10. Both enhancers and RNA data must have the same genome assembly

# Optional and arbitrary list with information and metadata that is stored
# within the GRaNIE object
objectMetadata.l = list(name = paste0("Macrophages_infected_primed"), file_peaks = file_peaks,
                        file_rna = file_RNA, file_sampleMetadata = file_sampleMetadata, genomeAssembly = genomeAssembly)

dir_output = "."

GRN = initializeGRN(objectMetadata = objectMetadata.l, outputFolder = dir_output,
                    genomeAssembly = genomeAssembly)

GRN

GRN = addData(GRN, counts_peaks = countsPeaks.df, normalization_peaks = "DESeq2_sizeFactors",
              idColumn_peaks = idColumn_peaks, counts_rna = countsRNA.df, normalization_rna = "limma_quantile",
              idColumn_RNA = idColumn_RNA, sampleMetadata = sampleMetadata.df, forceRerun = TRUE)

GRN

GRN = plotPCA_all(GRN, data = c("rna"), topn = 500, type = "normalized", plotAsPDF = FALSE,
                  pages = c(2, 3, 14), forceRerun = TRUE)

folder_TFBS_6TFs = "https://www.embl.de/download/zaugg/GRaNIE/TFBS_selected.zip"
# Download the zip of all TFBS files. Takes too long here, not executed
# therefore

download.file(folder_TFBS_6TFs, file.path("TFBS_selected.zip"), quiet = FALSE)

unzip(file.path("TFBS_selected.zip"), overwrite = TRUE)

motifFolder = tools::file_path_as_absolute("TFBS_selected")

GRN = addTFBS(GRN, motifFolder = motifFolder, TFs = "all", filesTFBSPattern = "_TFBS",
              fileEnding = ".bed.gz", forceRerun = TRUE)

GRN = overlapPeaksAndTFBS(GRN, nCores = 1, forceRerun = TRUE)
