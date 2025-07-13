# ICL-DESeq2_script
This is a DESeq2 script to run DESeq2 in Bash without interactively using R (or R Studio). Specifically, it is suited for running on HPC.

# TL;DR (Too long; Didn't read) version: #
**First time:**
```
conda create --name <env> --file [file path]/deseq2-list.txt
```
Make new directory for each analysis.
```
mkdir rnaseq_analysis_output
cp [file path to deseq2_script.R]/deseq2_script.R [file path]/rnaseq_analysis_output
```
Check file path for annotation folder and files (line 16).
Make a script file containing the following, with the right conda environment:
```
Rscript deseq2_script.R sample_table.txt
```
sample_table.txt example:

|SampleName|Treatment|filepath|
|----------|---------|--------|
|sample_DMSO_1|DMSO|[file path to salmon folder]/sample_DMSO_1|
|sample_DMSO_2|DMSO|[file path to salmon folder]/sample_DMSO_2|
|sample_ConditionA_1|ConditionA|[file path to salmon folder]/sample_ConditionA_1|
|sample_ConditionA_2|ConditionA|[file path to salmon folder]/sample_ConditionA_2|

File must be tab-separated.

# Main information
If this is your first time using this Rscript, please create a new conda environment before running the Rscript.\
The conda environment contains all the required R and bioconductor packages that are necessary for the Rscript to work.\
In bash, run the following line:
```
# Replace <env> with whatever name that you would like to call this conda environment, i.e. deseq2
# Replace [file.path] with wherever the deseq2-list.txt file was stored, preferably in the home directory, or in a directory containing all explicit lists of conda environments.
conda create --name <env> --file [file path]/deseq2-list.txt
```
Prior to running the DESeq2 script, [Salmon](https://github.com/COMBINE-lab/salmon) must be performed on the samples.\
For the script to work, each sample should be within its own directory, with the Salmon quant files stored within a subdirectory.\
If Salmon script was used, for each sample directory, a subdirectory called: new_salmon_quant/quant.sf (see line 71 of script).\
Additionally, the DESeq2 script only works with at least duplicates for each treatment.\
Due to the nature of how DESeq2 was written, this analysis only works for comparing between two treatment, i.e. DMSO vs Ponatinib, but not with 3 different conditions.\

**IMPORTANT NOTE!**\
**Always put DMSO (or vehicle) in the first row, otherwise, the script will always assume that the first row is the vehicle!**\
**This order is important, as DESeq2 will always compare the treatment to the vehicle/control/DMSO.**\
**This means that the fold change reported is compared to the vehicle/control/DMSO.**\

Please have the sample table in the following order for column: SampleName	Treatment	filepath\
Ensure that the header is spelt as shown above. Small typos can make the script fail!\
Additionally, the table should be tab-separated.\

Before running the Rscript, please check the file path to where the gene annotation files are stored.\
To check the file path, either open the script in R, or use nano in Bash and check line 16.\
It should say:
```
annotation_folder <- file.path("[file path to DESeq2 script]/DESeq2_script")
```

At the time of writing this script, the version of hg38 used is version 43, which was obtained from Gencode Genes.\
Go to https://www.gencodegenes.org/ and find the version of hg38 to download the comprehensive gene annotation GTF file.\
Should you wish to update to a newer gene annotation, please run the following in R and save the output and replace names of the files in lines 18, 25, 33 and 34.\
Replace the directory path in the line below with wherever you have the gtf file stored and where you would like save the sqlite file into.\
```
folder <- file.path("[file path to folder]")
hs_gft <- file.path(folder, "gencode.v43.annotation.gtf")
# Replace "Gencode_v43" with whatever version you have downloaded. #
txdb <- makeTxDbFromGFF(hs_gft, dataSource = "Gencode_v43", organism = "Homo sapiens")
txdb.filename <- gsub("gtf", "sqlite", hs_gft)
saveDb(txdb, txdb.filename)
# These four lines will generate a transcript to gene csv file. Similar to the sqlite file, the file will be saved in the same directory as written in the R object: folder
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
tx2gene.filename <- gsub("gtf", "tx2gene.csv", hs_gft)
write_csv(tx2gene, tx2gene.filename)
# The next five lines will generate a ENSEMBL gene annotation file for all genes and save it as an RDS file. As before, the file is saved into the same directory as "folder"
k2 <- keys(txdb, keytype = "GENEID")
k2 <- gsub("\\..*", "", k2)
gene_list <- select(org.Hs.eg.db, keys = k2, columns = c("SYMBOL","GENENAME"), keytype = "ENSEMBL")
gene_list <- gene_list[!duplicated(gene_list$ENSEMBL),]
saveRDS <- saveRDS(gene_list, file = file.path(folder, "RNAseq_annotated_gene_list.RDS"))
```
### Note for DESeq2 analysis
In the DESeq2 analysis, genes with very low count (of 5 or less) have been excluded from the analysis. Five is a conservative filtering parameter. The filtering value can be set to be more stringent.\
One reason is that the analysis is faster.\
Most likely, these genes are expressed at very low levels (or not at all) and outlier of expression from one replicate can skew the statistical analysis of what is truly a differentially expressed gene.\
i.e. Genes that are kept:\
|Gene|DMSO|DMSO|ConditionA|ConditionA|Total count|
|----|----|----|----------|----------|-----------|
|GeneA|0|0|3|3|6|
|GeneB|1|0|2|2|5|

i.e. Genes that are excluded:\
|Gene|DMSO|DMSO|ConditionA|ConditionA|Total count|
|----|----|----|----------|----------|-----------|
|GeneX|0|0|2|2|4|
|GeneY|0|0|0|0|0|

For more information on filtering, please read: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pre-filtering\
As such, in the final output file containing the normalised count for all genes, the genes with very low count of 5 or less have been excluded.\
The value of 5 is set for duplicate of each treatment. This value can be increased for more replicates.\
Additionally, an additional filtering can be made (in line 109) to keep total count by treatment group:
```
replicatesize <- X # replace X with the size of the replicate. #
keep5 <- rowSums(counts(dds) >= 5) >= replicatesize
```

One of the output file from DESeq2 analysis is a summary file: deseq2_summary_output.txt\
This txt summary file summarises all genes that passed the filtering threshold for a p-adjusted value of 0.05.\
A differentially expressed gene in the analysis is defined as log2 fold change of < -1 , or > 1, with a p-adjusted value < 0.05.

### Tips for changing aesthetic for plots generated with ggplot
All figures generated from DESeq2 are in A4 landscape orientation and with a dpi of 1200.\
Should the figures be changed for aesthetic reasons (or whatever the reason), open R studio, import the respective RDS file, and copy the script for the respective plot.\
For example, should the column label of the heatmap be changed to a different name that the default setting, run the following:
```
heatmap_data <- readRDS("[file path to DESeq2 output folder]/heatmap_data.RDS")
```
Add the following to change the label in pheatmap:
```
pheatmap(heatmap_data,show_rownames = FALSE, scale = "row", angle_col = 0, treeheight_col = 0, treeheight_row = 0, labels_col = c("DMSO\n(1)", "DMSO\n(2)", "ConditionA\n(1)", "ConditionA\n(2)")
```
Change the name with whatever is with each quotation.\
The \n syntax means that the anything after '\n' is printed on a new line.

Thank you for reading the document. Good luck with your RNA-seq analysis!
Last updated: 31/12/2024

# Output files from DESeq2 analysis
|File name|Description|
|---------|-----------|
|all_filtered_genes_DESeq2_stats.csv|This file contains the statistics (log2FoldChange & p-adj) for all genes after filtering.|
|dds.RDS|RDS file containing the DESeqDataSet of all genes (no filtering).|
|dds_filtered.RDS|RDS file containing the DESeqDataSet of genes after filtering (see above note).|
|ddsRes.RDS|RDS file containing the DESeqDataSet results, which has all statistics after DESeq2 analysis.|
|DEG_list_stats.csv|This file contains all DEG and the statistics (log2FoldChange & p-adj).|
|deseq2_summary_output.txt|This file contains an overall summary of ddsRes (with the p-adjusted value set to 0.05).|
|filtered_RNA_normalised_count.csv|This file contains the normalised count of all genes (after filtering).|
|heatmap_data.RDS|RDS file for the heatmap. Can be used to re-plot heatmap.|
|Heatmap_plot.pdf|PDF file of heatmap.|
|Heatmap_plot.png|PNG file of heatmap.|
|PCAplot.pdf|PDF file of PCA plot.|
|PCAplot.png|PNG file of PCA plot.|
|PCAplot_data.RDS|RDS file for the PCA plot. Can be used to re-plot PCA plot.|
|volcano_data.RDS|RDS file for the volcano plot. Can be used to re-plot volcano plot.|
|Volcano_plot.pdf|PDF file for volcano plot.|
|Volcano_plot.png|PNG file for volcano plot.|
