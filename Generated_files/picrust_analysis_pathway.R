##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Install packages if not already installed
for (pkg in pkgs) {
  print(pkg)
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# When prompted about updates, type "n" for none.

# After installing dependencies, install ggpicrust2
devtools::install_github("cafferychen777/ggpicrust2")

#### Load packages ####
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)

#### Import files and preparing tables ####
# Import the pathway PICrust2 abundance file
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read.delim(abundance_file, header = TRUE, skip = 1, sep = "\t")
abundance_data <- as.data.frame(abundance_data)

# Specify desktop path for figure export
desktop_path <- file.path(Sys.getenv("HOME"), "Desktop")

# Import your metadata file
metadata <- read_delim("ulcers_metadata.tsv")

# Remove NAs for your column of interest (subject)
metadata <- metadata[!is.na(metadata$subject), ]

# Filter metadata to only include samples where "Spaceflight" is either "Space Flight" or "Ground Control"
metadata <- metadata[metadata$Spaceflight %in% c("Space Flight", "Ground Control"), ]

# Set factor levels so "Ground Control" is the reference group
metadata$Spaceflight <- factor(metadata$Spaceflight, levels = c("Ground Control", "Space Flight"))

# Filtering abundance_data
# Use the metadata column "sample-id" for sample names
sample_names <- as.character(metadata$`sample-id`)
# Remove any occurrence of "X.OTU.ID" if it exists in sample_names (it shouldn’t, but just in case)
sample_names <- sample_names[!sample_names %in% "X.OTU.ID"]

# When subsetting the abundance table, keep the identifier column ("X.OTU.ID") plus the sample columns
abundance_data_filtered <- abundance_data[, colnames(abundance_data) %in% c("X.OTU.ID", sample_names)]

# (Optional) Inspect columns of abundance_data
print(colnames(abundance_data))

# Remove columns (samples) that have no data (zero counts) across features
abundance_data_filtered <- abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

# Extract sample names from abundance_data_filtered (excluding the identifier column "X.OTU.ID")
abun_samples <- colnames(abundance_data_filtered)[colnames(abundance_data_filtered) != "X.OTU.ID"]

# Ensure metadata only contains samples that are present in the abundance table
metadata <- metadata[metadata$`sample-id` %in% abun_samples, ]

#### DESeq Analysis ####
# Perform pathway DAA using DESeq2 method
abundance_daa_results_df <- pathway_daa(
  abundance = abundance_data_filtered %>% column_to_rownames("X.OTU.ID"),
  metadata = metadata,
  group = "subject",
  daa_method = "DESeq2"
)

# Annotate MetaCyc pathway for better descriptions
metacyc_daa_annotated_results_df <- pathway_annotation(
  pathway = "MetaCyc", 
  daa_results_df = abundance_daa_results_df, 
  ko_to_kegg = FALSE
)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_values < 0.05)

# Change the pathway column to description for the results
feature_desc <- inner_join(feature_with_p_0.05, metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature <- feature_desc$description
feature_desc <- feature_desc[, c(1:7)]
colnames(feature_desc) <- colnames(feature_with_p_0.05)

# Change the pathway column to description for the abundance table
abundance <- abundance_data_filtered %>% filter(`X.OTU.ID` %in% feature_with_p_0.05$feature)
colnames(abundance)[1] <- "feature"
abundance_desc <- inner_join(abundance, metacyc_daa_annotated_results_df, by = "feature")

# Overwrite the feature column with descriptive names (from "description")
abundance_desc$feature <- abundance_desc$description

# Instead of hard-coding removal of columns, subset to only the "feature" column and your sample columns
abundance_desc <- abundance_desc[, c("feature", sample_names)]
print(colnames(abundance_desc))

# -----------------------------
# Generate Figures
# -----------------------------

# Generate a heatmap
png(filename = file.path(desktop_path, "heatmap.png"), width = 1200, height = 600)
pathway_heatmap(
  abundance = abundance_desc %>% column_to_rownames("feature"),
  metadata = metadata,
  group = "Spaceflight"
)
dev.off()

# Add a new column "sample_name" to metadata (using the existing "sample-id" column)
metadata$sample_name <- metadata$`sample-id`

# Prepare abundance data for PCA by converting to rownames
abundance_for_pca <- abundance_data_filtered %>% column_to_rownames("X.OTU.ID")

# Remove any sample columns with zero variance (if any)
sample_sd <- apply(abundance_for_pca, 2, sd)
abundance_for_pca <- abundance_for_pca[, sample_sd > 0]

# Remove constant features (rows) that have zero variance across samples
feature_sd <- apply(abundance_for_pca, 1, sd)
abundance_for_pca <- abundance_for_pca[feature_sd > 0, ]

# Update metadata to include only the remaining samples
remaining_samples <- colnames(abundance_for_pca)
metadata <- metadata[metadata$sample_name %in% remaining_samples, ]


# Generate Bar Plot of log2FC from custom DESeq2 function
# Load the DESeq2 function
source("DESeq2_function.R")

# Run the function.
# Change "subject" below to the column you want to use for grouping if needed.
res <- DEseq2_function(abundance_data_filtered, metadata, "Spaceflight", id_col = "X.OTU.ID")

# Verify DESeq2 results
cat("DESeq2 results summary:\n")
print(summary(res))
cat("Unique feature names from DESeq2 results:\n")
print(unique(res$feature))

# Add a feature column if not already present (extra precaution)
res$feature <- rownames(res)

# Load stringr for name trimming
library(stringr)

# Create a new trimmed feature name column to match annotation names.
res$feature_trim <- str_remove(res$feature, "\\.\\.\\..*$")
cat("Unique trimmed feature names:\n")
print(unique(res$feature_trim))

# Verify the annotation table's feature names:
cat("Unique feature names in annotation table:\n")
print(unique(metacyc_daa_annotated_results_df$feature))

# Merge with pathway annotation using the trimmed feature names
res_desc <- inner_join(res, metacyc_daa_annotated_results_df, by = c("feature_trim" = "feature"))

# Check if the join worked:
cat("Merged results (first 6 rows):\n")
print(head(res_desc))
cat("Merged results summary:\n")
print(summary(res_desc))

# Optional: Drop unnecessary columns (adjust indices as needed)
# (Here, dropping columns 8:13 from the merged result)
cat("Column names before dropping extra columns:\n")
print(colnames(res_desc))
res_desc <- res_desc[, -c(8:13)]
cat("Column names after dropping extra columns:\n")
print(colnames(res_desc))

# Filter to only include significant pathways (with relaxed threshold)
sig_res <- res_desc %>% filter(pvalue < 0.05)
sig_res <- sig_res[order(sig_res$log2FoldChange), ]
nrow(sig_res)

# If sig_res is empty, revert to plotting all features (or adjust threshold further)
if(nrow(sig_res) == 0) {
  message("No significant pathways found at pvalue < 0.1, plotting all results.")
  sig_res <- res_desc
}

# Subset to the top 15 pathways based on pvalue (lowest p-value)
top15 <- sig_res %>% arrange(pvalue) %>% head(15)
# After subsetting:
# Check for NAs or duplicates in critical columns:
top15 %>% filter(is.na(description) | is.na(log2FoldChange))

# Generate the bar plot using the top 15 pathways (pathways with duplicate descriptions will not show twice)
bar_plot <- ggplot(data = top15, 
                   aes(y = reorder(description, log2FoldChange), x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Log2FoldChange", y = "Pathways")

# Optionally print the plot to the console to verify it visually:
print(bar_plot)

# Save the plot to file
png(filename = file.path(desktop_path, "bar_plot.png"), width = 800, height = 600)
print(bar_plot)
dev.off()


pca_result <- prcomp(t(abundance_for_pca), center = TRUE, scale. = TRUE)
summary(pca_result)
head(pca_result$x)

# Generate pathway PCA plot with the filtered data
png(filename = file.path(desktop_path, "pca_plot.png"), width = 800, height = 600)
pathway_pca(
  abundance = abundance_for_pca,
  metadata = metadata,
  group = "Spaceflight"
)
dev.off()

# Generate volcano plot for all pathways
# Remove rows with NA in pvalue or log2FoldChange
res_desc_clean <- res_desc %>% filter(!is.na(pvalue) & !is.na(log2FoldChange))
# Create a volcano plot using res_desc, which contains the DESeq2 results merged with pathway annotations.
volcano_plot <- ggplot(res_desc_clean, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = pvalue < 0.05)) +
  scale_color_manual(values = c("black", "red"), 
                     labels = c("Not Significant", "p < 0.05"),
                     name = "Significance") +
  labs(
    x = "log2(Fold Change) (Spaceflight vs. Ground Control)",
    y = "-log10(p-value)",
    title = "Volcano Plot: Differential Abundance of Pathways"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 12),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )


# Print volcano plot to console to verify it
print(volcano_plot)

# Save the volcano plot to a file on your desktop
png(filename = file.path(desktop_path, "volcano_plot.png"), width = 800, height = 600)
print(volcano_plot)
dev.off()

# Filter out rows with NA in pvalue and sort by pvalue
top_pathways <- res_desc %>% 
  filter(!is.na(pvalue)) %>% 
  arrange(pvalue) %>% 
  head(2)  # This returns the two most significant pathways

# Print the top two pathways
print(top_pathways)
# 1 palmitate biosynthesis II (bacteria and plants)
# 2    purine nucleobases degradation I (anaerobic)

# see how many pathways in total are changed
cat("Total significantly changed pathways (p < 0.05):", nrow(sig_res), "\n")
# Total significantly changed pathways (p < 0.05): 51
