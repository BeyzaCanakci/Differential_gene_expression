# Install necessary packages (if not already installed)
# install.packages(c("pheatmap", "viridis", "tidyverse"))
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("DESeq2", "EnhancedVolcano", "org.Hs.eg.db", "clusterProfiler", "pathview"))

# Load libraries
library(DESeq2)
library(tidyverse)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)
library(viridis)
library(clusterProfiler)
library(pathview)

# Read in count data
counts <- read.csv("GSE279527_read.counts.csv", sep = ",", header = TRUE, row.names = 1)

# Select relevant columns (first 3: control, last 3: treated)
counts <- select(counts, 2:7)

# Create metadata (sample information)
colData <- data.frame(condition = c(rep("control", 3), rep("treated", 3)))
rownames(colData) <- colnames(counts)

# Check if sample names match
stopifnot(all(rownames(colData) == colnames(counts)))

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~condition)

# Filter out lowly expressed genes
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]

# Set control as the reference level
dds$condition <- relevel(dds$condition, ref = "control")

# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)
res <- as.data.frame(res)

# MA Plot
png("plotMA.png", width = 900)
plotMA(dds, alpha = 0.01)  # alpha defines significance threshold
dev.off()

# PCA plot after variance stabilizing transformation
vstdata <- vst(dds, blind = FALSE)
plotPCA(vstdata, intgroup = "condition")

# Annotate results with gene SYMBOL and ENTREZ IDs
res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     keytype = "ENSEMBL",
                     column = "ENTREZID")
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     keytype = "ENSEMBL",
                     column = "SYMBOL")

# Volcano plot
EnhancedVolcano(res,
                x = "log2FoldChange",
                y = "padj",
                lab = res$symbol)

# Normalize counts
nor <- counts(dds, normalized = TRUE)

# Merge normalized counts with DESeq results
nor.res <- merge(res, nor, by = 0)

# Filter significant genes (padj < 0.05 and |log2FC| > 2)
significant.nor <- nor.res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::select(9:15)

# Set rownames for heatmap
rownames(significant.nor) <- make.names(significant.nor$symbol, unique = TRUE)
significant.nor$symbol <- NULL

# Draw heatmap for significant genes
pheatmap(log2(significant.nor + 1),  # log2 transformation for scale
         show_rownames = FALSE,
         color = viridis(50))

# Top 40 genes based on significance and strong fold change
top.genes <- nor.res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 3) %>%
  drop_na(symbol) %>%
  arrange(padj) %>%
  dplyr::select(9:15) %>%
  head(40)

rownames(top.genes) <- top.genes$symbol
top.genes$symbol <- NULL

# Heatmap of top 40 genes
pheatmap(log2(top.genes))

# Gene Ontology Enrichment
go.kegg <- nor.res %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2)

geneList <- go.kegg$Row.names

ego.mf <- enrichGO(gene = geneList,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "MF",
                   readable = TRUE)
barplot(ego.mf)

ego.cc <- enrichGO(gene = geneList,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "CC",
                   readable = TRUE)
barplot(ego.cc)

ego.bp <- enrichGO(gene = geneList,
                   OrgDb = "org.Hs.eg.db",
                   keyType = "ENSEMBL",
                   ont = "BP",
                   readable = TRUE)
barplot(ego.bp)

# KEGG Pathway Enrichment and Visualization
geneList.kegg <- go.kegg %>%
  drop_na(entrez) %>%
  select(entrez, log2FoldChange) %>%
  deframe()

gene <- names(geneList.kegg)

ekegg <- enrichKEGG(gene = gene,
                    organism = "hsa")  # Homo sapiens

ekegg <- as.data.frame(ekegg)

# Browse and visualize a specific KEGG pathway
browseKEGG(ekegg, "hsa04360")

pathview(gene.data = geneList.kegg,
         species = "hsa",
         pathway.id = "hsa04360",
         limit = list(gene = max(abs(geneList.kegg)), cpd = 1))
