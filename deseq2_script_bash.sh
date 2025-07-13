#!/bin/bash

# Load module containing the required R environment
# conda environment containing all R packages, i.e. tidyverse
module load r412

Rscript deseq2_script.R sample_table.txt
