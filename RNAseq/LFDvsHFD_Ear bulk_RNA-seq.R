renv::init()
renv::activate()
renv::status()

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
setwd("~/Downloads")
#prepare a metadata to include for information: 2 timepoints (D0 and D6) and 2 different diets (HFD and LFD)
sampleInfo <- data.frame(
  sampleName = c(
    paste0("D0-H", 1:5),
    paste0("D0-L", 1:5),
    paste0("D6-H", 1:5),
    paste0("D6-L", 1:5)
  ),
  time = rep(c("D0", "D0", "D6", "D6"), each = 5),
  diet = rep(c("HFD", "LFD", "HFD", "LFD"), each = 5)
)

#List all count files in the directory
count_files <- list.files(path = "/Users/scc/sophiecurio@googlemail.com - Google Drive/My Drive/Curio Science/BaparLab_Sequencing/GEO datasets/RNAseq", pattern = "_counts\\.txt$", full.names = TRUE)

# Initialize an empty list to store individual count data frames
count_list <- list()

# Read each count file and store the 7th column (read counts) in the list
for (i in seq_along(count_files)) {
  temp_data <- read.table(count_files[i], header=TRUE)
  count_list[[i]] <- temp_data[, 7, drop=FALSE]
  names(count_list)[i] <- basename(count_files[i])
}

# Combine all individual count data into a single count matrix
count_matrix <- do.call(cbind, count_list)

# Make sure to set row names as gene IDs from one of the individual count data
row.names(count_matrix) <- temp_data[, 1]


#-----------------------------compare LFD vs HFD at D0-----------------------------------------------------------------------
d0_samples <- sampleInfo$time == "D0"
count_matrix_d0 <- count_matrix[, d0_samples]
sampleInfo_d0 <- subset(sampleInfo, time == "D0")

#Create DESeq2 Object: Now, create the DESeq2 object using this subsetted count matrix and sample information. The design formula will only include the diet, since you're focusing on D0.

dds_D0 <- DESeqDataSetFromMatrix(
  countData = count_matrix_d0,
  colData = sampleInfo_d0,
  design = ~ diet
)

#Run DESeq2 Analysis: Finally, run DESeq2 on this subsetted object.
dds_D0 <- DESeq(dds_D0)

#-----------------
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("openxlsx")

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)

setwd("~/Downloads")

# Prepare metadata
sampleInfo <- data.frame(
  sampleName = c(
    paste0("D0-H", 1:5),
    paste0("D0-L", 1:5),
    paste0("D6-H", 1:5),
    paste0("D6-L", 1:5)
  ),
  time = rep(c("D0", "D0", "D6", "D6"), each = 5),
  diet = rep(c("HFD", "LFD", "HFD", "LFD"), each = 5)
)

# List all count files in the directory
count_files <- list.files(path = "/Users/scc/sophiecurio@googlemail.com - Google Drive/My Drive/Curio Science/BaparLab_Sequencing/GEO datasets/RNAseq", pattern = "_counts\\.txt$", full.names = TRUE)

# Initialize an empty list to store individual count data frames
count_list <- list()

# Read each count file and store the 7th column (read counts) in the list
for (i in seq_along(count_files)) {
  temp_data <- read.table(count_files[i], header = TRUE)
  count_list[[i]] <- temp_data[, 7, drop = FALSE]
  names(count_list)[i] <- basename(count_files[i])
}

# Combine all individual count data into a single count matrix
count_matrix <- do.call(cbind, count_list)

# Set row names as gene IDs from one of the individual count data
row.names(count_matrix) <- temp_data[, 1]

# Subset for D0 samples
d0_samples <- sampleInfo$time == "D0"
count_matrix_d0 <- count_matrix[, d0_samples]
sampleInfo_d0 <- subset(sampleInfo, time == "D0")

# Create DESeq2 object
dds_D0 <- DESeqDataSetFromMatrix(
  countData = count_matrix_d0,
  colData = sampleInfo_d0,
  design = ~ diet
)

# Run DESeq2 analysis
dds_D0 <- DESeq(dds_D0)

# Extract results for LFD vs HFD at D0
res_d0_diet <- results(dds_D0, contrast = c("diet", "HFD", "LFD"))

# Convert results to data frame
res_d0_diet_df <- as.data.frame(res_d0_diet)
res_d0_diet_df$gene_id <- rownames(res_d0_diet_df)

# Convert gene IDs to gene names
gene_names <- bitr(res_d0_diet_df$gene_id, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = org.Mm.eg.db)

# Handle 1:many mapping by keeping the first mapping
gene_names <- gene_names[!duplicated(gene_names$ENSEMBL), ]

# Merge gene names with results
res_d0_diet_df <- merge(res_d0_diet_df, gene_names, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)

# Order the data frame by adjusted p-value
res_d0_diet_df <- res_d0_diet_df[order(res_d0_diet_df$padj), ]

# Save results to Excel
write.xlsx(res_d0_diet_df, file = "DESeq2_LFD_vs_HFD_D0.xlsx", sheetName = "LFD_vs_HFD_D0", rowNames = FALSE)

save(res_d0_diet, file = "res_d0_diet.RData")


#---------------

#Extract Results: To get the differential gene expression results for D0 between HFD and LFD:
res_d0_diet <- results(dds_D0, contrast=c("diet", "HFD", "LFD"))

####visualization: Volcano plots##------------------------------------
mapping_file_path = "mart_export.csv"  ###need to convert ensembl ID to actual gene names
gene_mapping <- read.csv(mapping_file_path)
head(gene_mapping)

volcano_data <- as.data.frame(res_d0_diet)
volcano_data$Ensembl_ID <- rownames(volcano_data) 
volcano_data <- merge(volcano_data, gene_mapping, by="Ensembl_ID")
head(volcano_data)


# Define significance and fold-change thresholds
lpha <- 0.0001
log2FC_threshold <- 1

# Define colors based on significance and fold-change direction
volcano_data$color <- "grey"  # default color
volcano_data$color[volcano_data$padj < alpha & volcano_data$log2FoldChange > log2FC_threshold] <- "red"
volcano_data$color[volcano_data$padj < alpha & volcano_data$log2FoldChange < -log2FC_threshold] <- "blue"

# Create the volcano plot
library(ggplot2)
library(ggrepel)

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point() +
  geom_text_repel(data = subset(volcano_data, padj < alpha & abs(log2FoldChange) > log2FC_threshold),
                  aes(label = gene_name),
                  box.padding = 0.5) +
  scale_color_identity() +
  scale_x_continuous(limits = c(-2, 2)) +
  geom_hline(yintercept = -log10(alpha), linetype = "dashed") +  # significance threshold
  geom_vline(xintercept = log2FC_threshold, linetype = "dashed") +  # FC=1 threshold
  geom_vline(xintercept = -log2FC_threshold, linetype = "dashed") +  # FC=-1 threshold
  ggtitle("Volcano Plot for HFD over LFD") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value")

####visualization: Heatmap
library(pheatmap)

# Normalizing the counts
norm_counts <- counts(dds_D0, normalized=TRUE)

# Selecting the genes that are differentially expressed
significant_genes <- rownames(subset(res_d0_diet, padj < 0.05 & abs(log2FoldChange) > 1))

# Drawing the heatmap
heatmap_data <- norm_counts[significant_genes, ]
annotation_data <- data.frame(
  Group = factor(c(rep("HFD", 5), rep("LFD", 5)), levels = c("HFD", "LFD"))
)
# Define the color palette
annotation_colors = list(
  Group = c("HFD" = "red", "LFD" = "blue")
)
# Double-check that names match
stopifnot(all(names(annotation_colors$Group) %in% levels(annotation_data$Group)))

#Then use the following code to map Ensembl IDs to gene names:
library(biomaRt)

# Connect to the Ensembl database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get the mapping from Ensembl ID to Gene Name
genes <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
               filters = "ensembl_gene_id",
               values = rownames(heatmap_data), 
               mart = mart)

# Change the row names in heatmap_data
rownames(heatmap_data) <- genes$external_gene_name[match(rownames(heatmap_data), genes$ensembl_gene_id)]


rownames(annotation_data) <- colnames(heatmap_data)
pheatmap(heatmap_data,
         scale = "row",
         annotation_col = annotation_data,
         annotation_colors = annotation_colors,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_colnames = FALSE)
#--------------------------------GO enrichment analysis--D0--------------------- 
# Install and load necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")  # For mouse genes

library(clusterProfiler)
library(org.Mm.eg.db)
library(DESeq2)
library(ggplot2)

# Differential expression results (assuming 'dds_D0' is already defined and analyzed)
res_d0_diet <- results(dds_D0, contrast = c("diet", "HFD", "LFD"))

# Filter significant genes (adjust thresholds as necessary)
alpha <- 0.0001
log2FC_threshold <- 1

# Create a data frame of significant genes
sig_genes <- as.data.frame(res_d0_diet)
sig_genes <- sig_genes[which(sig_genes$padj < alpha), ]

# Convert ENSEMBL IDs to ENTREZ IDs
entrez_ids <- bitr(rownames(sig_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Merge with log2 fold changes
sig_genes <- merge(sig_genes, entrez_ids, by.x = "row.names", by.y = "ENSEMBL")
colnames(sig_genes)[1] <- "ENSEMBL"

# Perform GO enrichment analysis
go_enrichment <- enrichGO(gene          = sig_genes$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = "ENTREZID",
                          ont           = "BP",  # Can be "BP", "MF", or "CC"
                          pvalueCutoff  = 0.001,
                          qvalueCutoff  = 0.2)

# Extract the results
go_results <- as.data.frame(go_enrichment)

# Calculate the average log2 fold change for each GO term
go_term_to_genes <- bitr(unique(go_results$ID), fromType = "GOALL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
go_term_to_genes <- merge(go_term_to_genes, sig_genes, by.x = "ENTREZID", by.y = "ENTREZID")
go_avg_log2FC <- aggregate(log2FoldChange ~ GOALL, data = go_term_to_genes, FUN = mean)
go_results <- merge(go_results, go_avg_log2FC, by.x = "ID", by.y = "GOALL")

# Visualize the results
ggplot(go_results, aes(x = reorder(Description, log2FoldChange), y = log2FoldChange, fill = log2FoldChange)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  xlab("GO Term") +
  ylab("Average Log2 Fold Change") +
  ggtitle("GO Pathway Analysis: HFD vs LFD") +
  theme_minimal()

#-----------------------------compare LFD vs HFD at D6!!!!!-----------------------------------------------------------------------
d6_samples <- sampleInfo$time == "D6"
count_matrix_d6 <- count_matrix[, d6_samples]
sampleInfo_d6 <- subset(sampleInfo, time == "D6")

#Create DESeq2 Object: Now, create the DESeq2 object using this subsetted count matrix and sample information. The design formula will only include the diet, since you're focusing on D0.

dds_D6 <- DESeqDataSetFromMatrix(
  countData = count_matrix_d6,
  colData = sampleInfo_d6,
  design = ~ diet
)

#Run DESeq2 Analysis: Finally, run DESeq2 on this subsetted object.
dds_D6 <- DESeq(dds_D6)

#Extract Results: To get the differential gene expression results for D0 between HFD and LFD:
res_d6_diet <- results(dds_D6, contrast=c("diet", "HFD", "LFD"))


####visualization: Volcano plots##------------------------------------


volcano_data_D6 <- as.data.frame(res_d6_diet)
volcano_data_D6$Ensembl_ID <- rownames(volcano_data_D6)
head(volcano_data_D6$Ensembl_ID)

#Then use the following code to map Ensembl IDs to gene names:
library(biomaRt)

# Connect to the Ensembl database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get the mapping from Ensembl ID to Gene Name
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id",
                      values = rownames(volcano_data_D6), 
                      mart = mart)

volcano_data_D6 <- merge(volcano_data_D6, gene_mapping, by.x = "Ensembl_ID", by.y = "ensembl_gene_id")


head(volcano_data_D6)

library(ggplot2)
library(ggrepel)

# Define significance thresholds
alpha <- 0.0001
log2FC_threshold <- 1.5

# Create a new variable for the color gradient based on -log10(p-value)
volcano_data_D6$color_gradient <- ifelse(volcano_data_D6$padj >= alpha, -log10(volcano_data_D6$padj), NA)

# Base plot with gradient for non-significant genes
library(ggplot2)
library(ggrepel)


# Create the plot
p <- ggplot(volcano_data_D6, aes(x = log2FoldChange, y = -log10(padj))) +
  # Non-significant genes in a gray gradient
  geom_point(aes(color = log2FoldChange), size = 1.5, data = subset(volcano_data_D6, padj >= alpha)) +
  scale_color_gradient(low = "lightgrey", high = "darkgrey", na.value = NA) +
  # Significant genes with log2FoldChange > 1.5 in red
  geom_point(data = subset(volcano_data_D6, padj < alpha & log2FoldChange > log2FC_threshold),
             aes(x = log2FoldChange, y = -log10(padj)), color = "red", size = 2) +
  # Significant genes with log2FoldChange < -1.5 in blue
  geom_point(data = subset(volcano_data_D6, padj < alpha & log2FoldChange < -log2FC_threshold),
             aes(x = log2FoldChange, y = -log10(padj)), color = "blue", size = 2) +
  # Labels for significant genes
  geom_text_repel(data = subset(volcano_data_D6, padj < alpha & abs(log2FoldChange) > log2FC_threshold),
                  aes(label = external_gene_name), size = 5,
                  box.padding = unit(0.35, "lines")) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Volcano Plot for HFD over LFD", x = "Log2 Fold Change", y = "-Log10 P-value") +
  geom_hline(yintercept = -log10(alpha), linetype="dashed") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype="dashed")

# Print the plot
p




####visualization: Heatmap----------------------------------------------------
library(pheatmap)
library(DESeq2)
library(biomaRt)

# Normalizing the counts
norm_counts_D6 <- counts(dds_D6, normalized=TRUE)

# Selecting the genes that are differentially expressed
significant_genes_D6 <- rownames(subset(res_d6_diet, padj < 0.00001 & abs(log2FoldChange) > 1.5))

# Drawing the heatmap
heatmap_data_D6 <- norm_counts_D6[significant_genes_D6, ]
annotation_data_D6 <- data.frame(
  Group = factor(c(rep("HFD", 5), rep("LFD", 5)), levels = c("HFD", "LFD"))
)

# Define the color palette
annotation_colors_D6 = list(
  Group = c("LFD" = "blue", "HFD" = "red")
)

# Double-check that names match
stopifnot(all(names(annotation_colors_D6$Group) %in% levels(annotation_data_D6$Group)))

# Updated row names for annotation_data_D6
rownames(annotation_data_D6) <- colnames(heatmap_data_D6)

# Connect to the Ensembl database
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get the mapping from Ensembl ID to Gene Name
genes_D6 <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                  filters = "ensembl_gene_id",
                  values = rownames(heatmap_data_D6), 
                  mart = mart)

# Change the row names in heatmap_data_D6
rownames(heatmap_data_D6) <- genes_D6$external_gene_name[match(rownames(heatmap_data_D6), genes_D6$ensembl_gene_id)]

# Custom order of genes
custom_order = c("Serpina3b", "Serpinb6d", "Nxph3", "Pvalb",
                 "Il17a", "Il17f", "Osm", "Cxcl3", "Cxcl2", "Ccl3",  # Cytokines & Growth Factors
                 "Il1b", "Nlrp3", "Clec4e",  # Regulators & Modulators
                 "Vsig8", "Fpr1", "Trem1", "Trem3",  # Receptors
                 "Ptgs2", "Uox", "Cfap69",  # Enzymes
                 "AA467197", "Sprr2j-ps", "Klk14", "Slc26a4")  # Others

# Filter custom_order to include only genes present in heatmap_data_D6
custom_order <- custom_order[custom_order %in% rownames(heatmap_data_D6)]

# Get the indices for reordering
ordered_indices <- match(custom_order, rownames(heatmap_data_D6))

# Filter out NA values
ordered_indices <- ordered_indices[!is.na(ordered_indices)]

# Reorder the heatmap data
heatmap_data_D6_ordered <- heatmap_data_D6[ordered_indices, ]

# Update row names for annotation_data_D6
rownames(annotation_data_D6) <- colnames(heatmap_data_D6_ordered)

# Drawing the heatmap with increased label size
# [Your previous code for data preparation and reordering]

# Drawing the heatmap with clustering on rows disabled
pheatmap(heatmap_data_D6_ordered,
         scale = "row",
         annotation_col = annotation_data_D6,
         annotation_colors = annotation_colors_D6,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_rows = FALSE,  # Disable clustering on rows
         show_colnames = FALSE,
         fontsize = 12,  # Increased font size for column labels and legend
         fontsize_row = 20  # Increased font size for row labels
)



#--------------------------------GO enrichment analysis--D6--------------------- 
GO_genes_LFD_D6 <- subset(res_d6_diet, padj < 0.00001 & log2FoldChange < -1.5)

GO_genes_LFD_D6_ids <- rownames(GO_genes_LFD_D6)
# Perform GO enrichment analysis
GO_LFD_D6 <- enrichGO(gene         = GO_genes_LFD_D6_ids,
                      OrgDb        = org.Mm.eg.db, # replace with the OrgDb corresponding to your organism
                      keyType      = 'ENSEMBL',    # replace with the key type you are using
                      ont          = "ALL",        # options are "BP", "CC", "MF", and "ALL"
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

# Show the enrichment result
head(GO_LFD_D6)
result_df <- as.data.frame(GO_LFD_D6)
result_df_sorted <- result_df[order(result_df$pvalue), ]
GO_LFD_D6@result <- result_df_sorted

barplot(GO_LFD_D6, showCategory=8)

#------------HFD enrichment------
GO_genes_HFD_D6 <- subset(res_d6_diet, padj < 0.00001 & log2FoldChange > 1.5)

GO_genes_HFD_D6_ids <- rownames(GO_genes_HFD_D6)
# Perform GO enrichment analysis
GO_HFD_D6 <- enrichGO(gene         = GO_genes_HFD_D6_ids,
                      OrgDb        = org.Mm.eg.db, # replace with the OrgDb corresponding to your organism
                      keyType      = 'ENSEMBL',    # replace with the key type you are using
                      ont          = "BP",        # options are "BP", "CC", "MF", and "ALL"
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

# Show the enrichment result
head(GO_HFD_D6)
result_df <- as.data.frame(GO_HFD_D6)
result_df_sorted <- result_df[order(result_df$pvalue), ]
GO_HFD_D6@result <- result_df_sorted

# Take the top 20 most significant terms
result_df_top20 <- head(result_df_sorted, 20)

result_df_top20$significance_bin <- cut(-log10(result_df_top20$pvalue), breaks = 4, labels = c("Low", "Medium", "High", "Very High"))

library(RColorBrewer)
ggplot(result_df_top20, aes(x = reorder(Description, pvalue), y = -log10(pvalue), fill = significance_bin)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() +
  xlab("GO Term") +
  ylab("-log10(p-value)")


#------------looking at specific gene at those datasets 
#ENSMUSG00000000440 is PPARG in mouse
#ENSMUSG00000019577 is PDK4
#ENSMUSG00000025929 is Il17a
#ENSMUSG00000025888 is Casp1
#ENSMUSG00000032691 is Nlrp3
#ENSMUSG00000025746 is IL6
#ENSMUSG00000027398 is il1b
#Il1r1 ENSMUSG00000026072
#Gene: Il18 ENSMUSG00000039217

# Replace "PPARG" with the exact identifier used in your dataset
gene_id <- "ENSMUSG00000032691"

# Extract raw counts for PPARG at D0
counts_D0 <- counts(dds_D0)[gene_id, ]

# Extract raw counts for PPARG at D6
counts_D6 <- counts(dds_D6)[gene_id, ]

# Data frame for D0
plot_data_D0 <- data.frame(
  Diet = colData(dds_D0)$diet,
  Counts = counts_D0
)

# Data frame for D6
plot_data_D6 <- data.frame(
  Diet = colData(dds_D6)$diet,
  Counts = counts_D6
)

library(ggplot2)

# Plot for D0
ggplot(plot_data_D0, aes(x = Diet, y = Counts, color = Diet)) +
  geom_jitter(width = 0.2) +  # Adds some horizontal jitter to the points
  theme_minimal() +
  labs(title = "PPARG Raw Counts at D0", x = "Diet", y = "Raw Counts")

# Plot for D6
ggplot(plot_data_D6, aes(x = Diet, y = Counts, color = Diet)) +
  geom_jitter(width = 0.2) +  # Adds some horizontal jitter to the points
  theme_minimal() +
  labs(title = "PPARG Raw Counts at D6", x = "Diet", y = "Raw Counts")



