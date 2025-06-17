## Installation
#BiocManager::install("tximport")
library("tximport")
BiocManager::install("pasilla")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("DEGreport")
BiocManager::install("rtracklayer")
BiocManager::install("Glimma")
ainstall.packages("pheatmap")


## Initiate libraries.
library("pasilla")
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("pasilla")
library("DEGreport")
library("rtracklayer")
library("tximport")
library("Glimma")
library("edgeR")


# Create the sample_id vector
sample_ids <- c("st_21_22_1_", "st_21_22_2_", "st_21_22_3_",
                "st_22_23_1_", "st_22_23_2_", "st_22_23_3_",
                "st_24_25_1_", "st_24_25_2_", "st_24_25_3_",
                "st_25_25_1_", "st_25_25_2_", "st_25_25_3_",
                "st_27_28_1_", "st_27_28_2_", "st_27_28_3_",
                "st_29_29_1_", "st_29_29_2_", "st_29_29_3_",
                "st_30_31_1_", "st_30_31_2_", "st_30_31_3_",
                "st_31_32_1_", "st_31_32_2_", "st_31_32_3_",
                "st_32_33_1_", "st_32_33_2_", "st_32_33_3_",
                "st_33_34_1_", "st_33_34_2_", "st_33_34_3_")

# Create the samples data frame
samples <- data.frame(
  sample_id = sample_ids,
  stages = sapply(strsplit(sample_ids, "_"), function(x) paste(x[2], x[3], sep = "_"))
)


base_dir <- "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/salmon_quant"
files <- file.path(base_dir, paste0(samples$sample_id, "_quant"), "quant.sf")
names(files) <- samples$sample_id
# Convert stages to factor for DESeq2
samples$stages <- factor(samples$stages)



# Set row names (important for DESeq2)
rownames(samples) <- samples$sample_id
# Display the sample information
print(samples[,c("sample_id", "stages")])

gtf <- import("/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/GCF_040954835.1_fLepOcu1.hap2_genomic.gtf")

tx2gene <- read_csv("/mnt/research/FishEvoDevoGeno/raw_reads_4_tony/tx2gene.csv") 
clean_salmon_ids <- function(file_path) {
  # Read the salmon quant.sf file
  quant_data <- read.table(file_path, header = TRUE, sep = "\t")
  
  # Extract clean XR IDs
  quant_data$Name <- str_extract(quant_data$Name, "XR_\\d+\\.\\d+")
  
  # Remove rows where extraction failed
  quant_data <- quant_data[!is.na(quant_data$Name), ]
  
  # Write back to file (or create new file)
  write.table(quant_data, file = paste0(file_path, "_cleaned"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  return(paste0(file_path, "_cleaned"))
}
files_cleaned <- sapply(files, clean_salmon_ids)
txi <- tximport(files_cleaned, 
                type = "salmon", 
                tx2gene = tx2gene)
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ stages) 


## Differential expression analysis
dds <- DESeq(ddsTxi)
colnames(dds)
res <- results(dds)
resultsNames(dds)

dds_1 <- DESeq(ddsTxi, quiet=TRUE)
glimmaMDS(dds)
glimmaMA(dds_1,groups=dds_1$stages)
glimmaVolcano(dds_1,groups=dds_1$stages)
ncol(counts(dds_1))


######################################################################
## Log fold change shrinkage for visualization and ranking
######################################################################
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="stages_33_34_vs_21_22", type="apeglm")
resLFC

plotMA(resLFC, ylim=c(-2,2))


res_specific <- results(dds, contrast=c("stages", "29_29", "30_31"))
glimmaVolcano(dds_1, dge=res_specific, groups=dds_1$stages)

## ggplot2 (looking cute)
library(ggplot2)
library(ggrepel)
library(dplyr)

# Assuming you have your DESeq2 results in 'res' or 'res_shrunk'
# Replace 'res' with your actual results object name
res_df <- as.data.frame(res)  # or just 'res' if not using shrinkage

res_df$gene <- rownames(res_df)

# Remove rows with NA values
res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]

# Set significance thresholds
padj_threshold <- 0.05
lfc_threshold <- 1

# Create categories for coloring
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange > lfc_threshold] <- "Upregulated"
res_df$significance[res_df$padj < padj_threshold & res_df$log2FoldChange < -lfc_threshold] <- "Downregulated"
res_df$significance[res_df$padj < padj_threshold & abs(res_df$log2FoldChange) < lfc_threshold] <- "Significant"

# Define colors
colors <- c("Upregulated" = "#E31A1C", 
            "Downregulated" = "#1F78B4", 
            "Significant" = "#33A02C",
            "Not significant" = "grey60")

# Create the volcano plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1) +
  scale_color_manual(values = colors) +
  
  # Add threshold lines
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.7) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", alpha = 0.7) +
  
  # Customize the plot
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    legend.position = "none",  # Remove legend like in the original
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  
  labs(
    title = "RNAseq volcano plot",
    x = "Log 2 fold change",
    y = "-Log 10 P"
  ) +
  
  # Set axis limits (adjust based on your data)
  xlim(c(-6, 6)) +
  ylim(c(0, max(-log10(res_df$padj), na.rm = TRUE) + 1))

# Add gene labels for top significant genes
# Select top genes to label (adjust criteria as needed)
top_genes <- res_df %>%
  filter(padj < 0.001 & abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  head(20)  # Top 20 most significant genes

# If you want to label specific genes, you can also do:
# top_genes <- res_df %>%
#   filter(gene %in% c("GENE1", "GENE2", "GENE3", ...))

volcano_plot_labeled <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    box.padding = 0.3,
    point.padding = 0.3,
    max.overlaps = 20,
    segment.color = "black",
    segment.size = 0.2
  )


# Display the plot
print(volcano_plot_labeled)

## extracting transformed values
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
##############################################################
## heatmap of the sample-to-sample distances and PCA
##############################################################
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$stages, vsd$sample_id, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
colnames(colData(dds))

## Principal component plot of the samples
plotPCA(vsd, intgroup=c("stages"))

## Mean-Variance QC plots
counts<-counts(dds,normalized = TRUE)
design <- as.data.frame(colData(dds))
degQC(counts, design[["stages"]], pvalue = res[["pvalue"]])

## Covariates correlation with metrics
cor <- degCorCov(colData(dds))


###############################################################################
### pairwise comparison: 
###############################################################################
# Function to create volcano plot
create_volcano_plot <- function(res_df, comparison_name, 
                               padj_threshold = 0.05, 
                               lfc_threshold = 1) {
  
  # Define colors for significance
  colors <- c(
    "Up-regulated" = "red",
    "Down-regulated" = "blue", 
    "Not significant" = "grey"
  )
  
  # Add significance column
  res_df$significance <- ifelse(is.na(res_df$padj), "Not significant",
                               ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange > lfc_threshold, "Up-regulated",
                                     ifelse(res_df$padj < padj_threshold & res_df$log2FoldChange < -lfc_threshold, "Down-regulated", 
                                           "Not significant")))
  
  # Create the volcano plot
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = significance), alpha = 0.6, size = 1) +
    scale_color_manual(values = colors) +
    
    # Add threshold lines
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), linetype = "dashed", color = "black", alpha = 0.7) +
    
    # Customize the plot
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 16, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    
    labs(
      title = paste("Volcano Plot:", comparison_name),
      x = "Log 2 fold change",
      y = "-Log 10 adjusted P-value"
    ) +
    
    # Set axis limits
    xlim(c(-6, 6)) +
    ylim(c(0, max(-log10(res_df$padj), na.rm = TRUE) + 1))
  
  # Add gene labels for top significant genes
  top_genes <- res_df %>%
    filter(padj < 0.001 & abs(log2FoldChange) > 2) %>%
    arrange(padj) %>%
    head(20)
  
  if(nrow(top_genes) > 0) {
    volcano_plot <- volcano_plot +
      geom_text_repel(
        data = top_genes,
        aes(label = rownames(top_genes)),  # or use gene column if available
        size = 3,
        box.padding = 0.3,
        point.padding = 0.3,
        max.overlaps = 20,
        segment.color = "black",
        segment.size = 0.2
      )
  }
  
  return(volcano_plot)
}

# This compares every stage against every other stage
all_pairwise_comparisons <- function(dds) {
  stage_levels <- levels(dds$stages)
  comparison_results <- list()
  
  for(i in 1:length(stage_levels)) {
    for(j in 1:length(stage_levels)) {
      if(i != j) {  # Don't compare a stage to itself
        stage1 <- stage_levels[i]
        stage2 <- stage_levels[j]
        
        comparison_name <- paste(stage1, "vs", stage2)
        
        # Use contrast to specify the comparison
        res <- results(dds, contrast = c("stages", stage1, stage2))
        
        # Convert to data frame
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        
        # Remove rows with NA padj
        res_df <- res_df[!is.na(res_df$padj), ]
        
        # Create volcano plot
        volcano_plot <- create_volcano_plot(res_df, comparison_name)
        
        # Store results
        comparison_results[[comparison_name]] <- list(
          results = res_df,
          plot = volcano_plot
        )
      }
    }
  }
  
  return(comparison_results)
}

#For all pairwise comparisons (warning: many comparisons!)
all_results <- all_pairwise_comparisons(dds)


# Create the directory structure
base_path <- "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony"
plot_folder <- file.path(base_path, "plot", "RNA", "pairwise_comparison")

# Create directories if they don't exist
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Check if directory was created successfully
if(dir.exists(plot_folder)) {
  cat("Successfully created directory:", plot_folder, "\n")
} else {
  cat("Failed to create directory:", plot_folder, "\n")
  stop("Directory creation failed")
}

# Function to clean up comparison names for filenames
clean_filename <- function(comparison_name) {
  # Replace spaces and special characters with underscores
  filename <- gsub(" vs ", "_vs_", comparison_name)
  filename <- gsub("[^A-Za-z0-9_]", "_", filename)
  # Remove multiple consecutive underscores
  filename <- gsub("_{2,}", "_", filename)
  # Remove leading/trailing underscores
  filename <- gsub("^_|_$", "", filename)
  return(filename)
}

# Save all pairwise comparison plots
cat("Saving", length(all_results), "pairwise comparison plots...\n")

saved_files <- c()
failed_saves <- c()

for(comparison_name in names(all_results)) {
  tryCatch({
    # Clean the filename
    clean_name <- clean_filename(comparison_name)
    filename <- paste0("volcano_", clean_name, ".png")
    filepath <- file.path(plot_folder, filename)
    
    # Save the plot
    ggsave(filepath, 
           all_results[[comparison_name]]$plot, 
           width = 10, height = 8, dpi = 300)
    
    saved_files <- c(saved_files, filename)
    cat("✓ Saved:", filename, "\n")
    
  }, error = function(e) {
    failed_saves <- c(failed_saves, comparison_name)
    cat("✗ Failed to save:", comparison_name, "- Error:", e$message, "\n")
  })
}

# Create the directory structure
base_path <- "/mnt/research/FishEvoDevoGeno/raw_reads_4_tony"
plot_folder <- file.path(base_path, "plot", "RNA", "pairwise_comparison")

# Create directories if they don't exist
dir.create(plot_folder, recursive = TRUE, showWarnings = FALSE)

# Check if directory was created successfully
if(dir.exists(plot_folder)) {
  cat("Successfully created directory:", plot_folder, "\n")
} else {
  cat("Failed to create directory:", plot_folder, "\n")
  stop("Directory creation failed")
}

# Function to clean up comparison names for filenames
clean_filename <- function(comparison_name) {
  # Replace spaces and special characters with underscores
  filename <- gsub(" vs ", "_vs_", comparison_name)
  filename <- gsub("[^A-Za-z0-9_]", "_", filename)
  # Remove multiple consecutive underscores
  filename <- gsub("_{2,}", "_", filename)
  # Remove leading/trailing underscores
  filename <- gsub("^_|_$", "", filename)
  return(filename)
}

# Save all pairwise comparison plots
cat("Saving", length(all_results), "pairwise comparison plots...\n")

saved_files <- c()
failed_saves <- c()

for(comparison_name in names(all_results)) {
  tryCatch({
    # Clean the filename
    clean_name <- clean_filename(comparison_name)
    filename <- paste0("volcano_", clean_name, ".png")
    filepath <- file.path(plot_folder, filename)
    
    # Save the plot
    ggsave(filepath, 
           all_results[[comparison_name]]$plot, 
           width = 10, height = 8, dpi = 300)
    
    saved_files <- c(saved_files, filename)
    cat("✓ Saved:", filename, "\n")
    
  }, error = function(e) {
    failed_saves <- c(failed_saves, comparison_name)
    cat("✗ Failed to save:", comparison_name, "- Error:", e$message, "\n")
  })
}

# To access a specific comparison:
# sequential_results[["st_22_23 vs st_21_22"]]$plot  # The volcano plot
# sequential_results[["st_22_23 vs st_21_22"]]$results  # The results data frame

# To save plots:
# ggsave("volcano_st_22_23_vs_st_21_22.png", 
#        sequential_results[["st_22_23 vs st_21_22"]]$plot, 
#        width = 10, height = 8, dpi = 300)

# Example of how to modify your existing code to work with specific comparisons:
# Replace 'res_df' in your existing volcano plot code with:
# res_df <- sequential_results[["st_22_23 vs st_21_22"]]$results


