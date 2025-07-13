#! /usr/bin/env Rscript

# Libraries required for RNA-seq analysis ----
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tximport)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(GenomicFeatures)))
suppressWarnings(suppressMessages(library(org.Hs.eg.db)))
suppressWarnings(suppressMessages(library(magrittr)))
suppressWarnings(suppressMessages(library(DESeq2)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(pheatmap)))
suppressWarnings(suppressMessages(library(biomaRt)))

# Gene annotation files (Change the file path to match where the annotation files are stored)
annotation_folder <- file.path("[file path to DESeq2]/DESeq2_script")
# txdb 
txdb.filename <- file.path(annotation_folder,"gencode.v43.annotation.sqlite")
txdb <- loadDb(txdb.filename)
if (!file.exists(txdb.filename)) {
  stop("The txdb file does not exist. Please check if the file exist in the file path: ",txdb.filename)
}
# Generate a tx2gene object for tximport
k <- keys(txdb, keytype = "TXNAME")
tx2gene.filename <- file.path(annotation_folder, "gencode.v43.annotation.tx2gene.csv")
tx2gene <- read_csv(tx2gene.filename)
if (!file.exists(tx2gene.filename)) {
  stop("The tx2gene file does not exist. Please check if the file exists in the file path: ", tx2gene.filename)
}
k2 <- keys(txdb, keytype = "GENEID")
k2 <- gsub("\\..*", "", k2)
# Gene list (a table of transcript and the corresponding gene)
gene_list_filepath <- file.path(annotation_folder,"RNAseq_annotated_gene_list.RDS")
gene_list <- readRDS(file.path(annotation_folder,"RNAseq_annotated_gene_list.RDS"))
if (!file.exists(gene_list_filepath)) {
  stop("The gene list file does not exist. Please check if the file exists in the file path: ", gene_list_filepath)
}

# Set the output directory for all files to be generated
output_dir <- getwd()

# Parse bash command line argument into R
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample table. Usage: Rscript deseq2_script.R <path_to_sample_table.txt>", call. = FALSE)
}

# Input file
input_file <- args[1]
cat("Input file:", input_file, "\n")

# Check that the file exists (i.e. in the right file path)
if (!file.exists(input_file)) {
  stop("The input file does not exist. Please check the file path. \n 
       Example: Rscript deseq2_script.R ./sample_file.txt")
}

# Read the input file. Input file in the command line should be in tab-delimited format!
# The input file should be a table containing the following columns:
#   - SampleName: Identifier for the sample (each sample should have a unique identifier)
#   - Treatment: Treatment of samples with either a vehicle, or treatment
#   - filepath: absolute file path to where the Salmon Quant files are stored
sample_table <- read.table(input_file, header = TRUE, sep = "\t")

# Validate the structure of sample_table
if (!all(c("SampleName","Treatment","filepath") %in% colnames(sample_table))) {
  stop("The input file must have the following columns: 'SampleName', 'Treatment', 'filepath'.")
}

# Import Salmon quant.sf files into R 
sample_files <- paste0(pull(sample_table,"filepath"),"/new_salmon_quant/quant.sf")
names(sample_files) <- pull(sample_table,"SampleName")

# Use tximport to convert quant.sf files into an R object
txi <- tximport(sample_files,
                type = "salmon",
                tx2gene = tx2gene)
if (!exists("txi") || is.null(txi)) {
  stop("Failed to create tximport object. Please check file path to quant.sf and tx2gene.")
}

# Setting factors and levels of treatment from the sample_table
# This step is important to ensure that the first factor is always the control group
# The vehicle sample is always set as first
sample_table <- as.data.frame(sample_table)
control_level <- as.character(sample_table$Treatment[1])
treatment_level <- unique(sample_table$Treatment[sample_table$Treatment != control_level])[1]
table_condition <- sample_table$Treatment
table_condition <- factor(table_condition, levels = c(control_level,treatment_level))
sample_table$Treatment <- table_condition
cat("The control group (first level) is set as:",control_level,"\n")
cat("The levels in the 'Treatment' columns are:\n")
cat(paste(levels(sample_table$Treatment), collapse = ", "), "\n")

# Convert tximport file into a DESeqDataSet
dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = sample_table,
                                design = ~Treatment)
if (!exists("dds") || is.null(dds)) {
  stop("Failed to create DESeqDataSet object. Please check tximport data.")
}

# save dds object back out of R
saveRDS(dds, file = file.path(output_dir,"dds.RDS"))
cat("The DESeqDataSet object (dds) has been saved to:", file.path(output_dir,"dds.RDS"), "\n")

# Filter out genes with very low count from the analysis
# See DESeq2 tutorial on pre-filtering for why it is normally recommended for DESeq2 analysis
keep5 <- rowSums(counts(dds)) >= 5
dds_filtered <- dds[keep5,]
dds_filtered <- DESeq(dds_filtered)
if (!exists("dds_filtered") || is.null(dds_filtered)) {
 stop("Failed to filter DESeqDataSet object. Check filtering criteria.") 
}
saveRDS(dds_filtered, file = file.path(output_dir,"dds_filtered.RDS"))
cat("The filtered DESeqDataSet object (dds_filtered) has been processed by DESeq2 and saved to:", file.path(output_dir,"dds_filtered.RDS"), "\n")

# Annotation of the normalised gene count after performing DESeq2
normalised_count_filtered <- as.data.frame(counts(dds_filtered, normalized = TRUE))
normalised_count_filtered_annotated <- rownames_to_column(normalised_count_filtered, var = "ENSEMBL_GENEID")
normalised_count_filtered_annotated$ENSEMBL_GENEID <- gsub("\\..*","",normalised_count_filtered_annotated$ENSEMBL_GENEID)
normalised_count_filtered_annotated <- dplyr::left_join(normalised_count_filtered_annotated, gene_list,
                                                        by = c("ENSEMBL_GENEID"="ENSEMBL"))
normalised_count_filtered_annotated <- normalised_count_filtered_annotated %>%
  dplyr::select("ENSEMBL_GENEID","SYMBOL",everything(),"GENENAME")
write.csv(normalised_count_filtered_annotated, row.names = FALSE,
          file.path(output_dir,"filtered_RNA_normalised_count.csv"))
cat("The normalised gene count (normalised_count_filtered_annotated) has been saved to:",file.path(output_dir,"filtered_RNA_normalised_count.csv"),"\n")

# Results of DESeq2 analysis for differentially expressed genes
ddsRes <- results(dds_filtered)
ddsRes005 <- results(dds_filtered, alpha = 0.05)
summary(ddsRes005)
capture.output(summary(ddsRes005), file = file.path(output_dir,"deseq2_summary_output.txt"))
cat("The summary output of DESeq2 results has been saved to:", file.path(output_dir,"deseq2_summary_output.txt"),"\n")
saveRDS(ddsRes, file = file.path(output_dir,"ddsRes.RDS"))
cat("The dds result object (ddsRes) has been saved to:", file.path(output_dir,"ddsRes.RDS"), "\n")

# Analysis of differentially expressed genes
# The following ddsRes_df has statistics on all genes from dds_filtered.
ddsRes_df <- as.data.frame(ddsRes)
all_ddsRes_annotated <- rownames_to_column(ddsRes_df, var = "ENSEMBL_GENEID")
all_ddsRes_annotated$ENSEMBL_GENEID <- gsub("\\..*","",all_ddsRes_annotated$ENSEMBL_GENEID)
all_ddsRes_annotated <- dplyr::left_join(all_ddsRes_annotated, gene_list,
                                         by = c("ENSEMBL_GENEID"="ENSEMBL"))
all_ddsRes_annotated <- all_ddsRes_annotated %>%
  dplyr::select("ENSEMBL_GENEID","SYMBOL",everything(),"GENENAME")
write.csv(all_ddsRes_annotated, row.names = FALSE,
          file.path(output_dir,"all_filtered_genes_DESeq2_stats.csv"))
cat("For all genes with a total count greater than 5 across all samples, the DESeq2 statistics (all_ddsRes_annotated) has been saved to:",
    file.path(output_dir,"all_filtered_genes_DESeq2_stats.csv"), "\n")

# The following filtering list out the differentially expressed genes (DEGs).
# In the following filtering, a DEG has a p-adjusted value < 0.05 and a log2 fold change greater than 1.
ddsRes_df_filter1 <- ddsRes_df[complete.cases(ddsRes_df),]
ddsRes_df_filter2 <- ddsRes_df_filter1[ddsRes_df_filter1$padj < 0.05,]
ddsRes_df_filter3 <- ddsRes_df_filter2[abs(ddsRes_df_filter2$log2FoldChange) > 1,]

DEG_list <- as.data.frame(ddsRes_df_filter3)
DEG_list_annotated <- rownames_to_column(DEG_list, var = "ENSEMBL_GENEID")
DEG_list_annotated$ENSEMBL_GENEID <- gsub("\\..*","",DEG_list_annotated$ENSEMBL_GENEID)
DEG_list_annotated <- dplyr::left_join(DEG_list_annotated, gene_list,
                                       by = c("ENSEMBL_GENEID"="ENSEMBL"))
DEG_list_annotated <- DEG_list_annotated %>%
  dplyr::select("ENSEMBL_GENEID","SYMBOL",everything(),"GENENAME")
write.csv(DEG_list_annotated, row.names = FALSE,
          file.path(output_dir,"DEG_list_stats.csv"))
cat("All DEGs (DEG_list_annotated) have been saved to:", file.path(output_dir,"DEG_list_stats.csv"), "\n")

# Plot PCA
dds_rld <- rlog(dds_filtered)
saveRDS(dds_rld, file = file.path(output_dir,"PCAplot_data.RDS"))
cat("Rlog transformation was performed on dds_filtered (dds_rld) and saved to:", file.path(output_dir,"PCAplot_data.RDS"), "\n")
cat("If there are changes to be made to the PCA plot, please import PCAplot_data.RDS into R to make modifications with ggplot.", "\n")

pca_plot <- plotPCA(dds_rld, intgroup = c("Treatment")) + theme_bw() +
  theme(plot.background = element_blank(), panel.grid = element_blank(), aspect.ratio = 1)
# PNG file format (A4 landscape) for PCA plot
pca_png <- file.path(output_dir,"PCAplot.png")
ggsave(filename = pca_png, plot = pca_plot, device = "png", width = 29.7, height = 21.0, units = "cm", dpi = 1200)
# PDF file format (A4 landscape) for PCA plot
pca_pdf <- file.path(output_dir,"PCAplot.pdf")
ggsave(filename = pca_pdf, plot = pca_plot, device = "pdf", width = 29.7, height = 21.0, units = "cm", dpi = 1200)

# Plot heatmap
DEG_genes <- rownames(DEG_list)
heatmap_data <- assay(dds_rld)[DEG_genes,]
saveRDS(heatmap_data, file = file.path(output_dir,"heatmap_data.RDS"))
cat("Heatmap data (heatmap_data) is saved to:", file.path(output_dir,"heatmap_data.RDS"), "\n")

heatmap_plot <- pheatmap(heatmap_data,show_rownames = FALSE, scale = "row", angle_col = 0, 
                         treeheight_col = 0, treeheight_row = 0, silent = TRUE)
# PNG file format (A4 landscape) for heatmap plot
heatmap_png <- file.path(output_dir,"Heatmap_plot.png")
ggsave(filename = heatmap_png, plot = heatmap_plot, device = "png", width = 29.7, height = 21.0, units = "cm", dpi = 1200)
# PDF file format (A4 landscape) for heatmap plot
heatmap_pdf <- file.path(output_dir,"Heatmap_plot.pdf")
ggsave(filename = heatmap_pdf, plot = heatmap_plot, device = "pdf", width = 29.7, height = 21.0, units = "cm", dpi = 1200)

# Plot volcano plot
volcano_df <- ddsRes_df_filter1
volcano_df <- volcano_df %>% dplyr::mutate(
  expression = case_when(log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
                         log2FoldChange > 1 & padj < 0.05 ~ "Up-regulated",
                         TRUE ~ "Unchanged"))
saveRDS(volcano_df, file = file.path(output_dir,"volcano_data.RDS"))
cat("Volcano plot data (volcano_df) is saved to: ", file.path(output_dir,"volcano_data.RDS"), "\n")

volcano_plot <- ggplot(volcano_df, aes(x=log2FoldChange, y=-log10(padj))) + theme_bw() +
  geom_point(size = 3, aes(colour = expression)) +
  scale_colour_manual(values = c("Unchanged"="black","Down-regulated"="red","Up-regulated"="green")) +
  geom_vline(xintercept = -1, linetype = "dotted") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  theme(legend.position = "none", panel.grid = element_blank())
# PNG file format (A4 landscape) for volcano plot
volcano_png <- file.path(output_dir,"Volcano_plot.png")
ggsave(filename = volcano_png, plot = volcano_plot, device = "png", width = 29.7, height = 21.0, units = "cm", dpi = 1200)
# PDF file format (A4 landscape) for volcano plot
volcano_pdf <- file.path(output_dir,"Volcano_plot.pdf")
ggsave(filename = volcano_pdf, plot = volcano_plot, device = "pdf", width = 29.7, height = 21.0, units = "cm", dpi = 1200)

cat("DESeq2 analysis completed! \n")
