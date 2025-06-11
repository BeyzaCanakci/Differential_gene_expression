
# Differential Gene Expression Analysis (Bulk RNA-Seq)

This repository contains code for performing differential gene expression analysis using bulk RNA-Seq data in R.

## Overview

The main script, differential_gene_expression_analysis_bulkRNA.R, includes all the necessary steps for processing RNA-Seq count data, performing normalization, and identifying differentially expressed genes between experimental conditions.

## Features

- Data import and preprocessing  
- Normalization of raw counts  
- Differential expression analysis using popular R packages (e.g., DESeq2, edgeR, limma)  
- Visualization of results (e.g., heatmaps, volcano plots)  
- Export of significant gene lists

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/BeyzaCanakci/Differential_gene_expression.git
   ```

2. Open the R script (`differential_gene_expression_analysis_bulkRNA.R`) in RStudio or your preferred R environment.

3. Update the file paths and parameters as needed for your data.

4. Run the script to perform the analysis.

## Requirements

- R (version 4.0 or higher recommended)
- R packages:
  - DESeq2
  - edgeR
  - limma
  - ggplot2
  - pheatmap
  - (and any other packages used in the script)

Install required packages in R:
```R
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "edgeR", "limma", "pheatmap"))
install.packages("ggplot2")
```

## Input Data

- Raw count matrix (genes x samples)
- Sample information/metadata (e.g., conditions, replicates)

## Output

- Tables of differentially expressed genes
- Plots for quality control and results visualization


## Contact

For questions or suggestions, please contact [Beyza Canakci](https://github.com/BeyzaCanakci).

---
