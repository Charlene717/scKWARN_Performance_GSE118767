## Quality control for Seurat Object

# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

if(!require("Seurat")) install.packages("Seurat"); library(Seurat)


seurat_list <- lapply(seurat_list, function(x) {
  # # Way1: Calculate the percentage of mitochondrial genes using the Seurat function
  # x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")

  # Way2: Manually calculate the percentage of mitochondrial genes
  total_counts_per_cell <- Matrix::colSums(x[["RNA"]]@counts)
  mito_genes <- grep("^MT-", rownames(x[["RNA"]]@counts), value = TRUE)
  x[["percent.mt"]] <- Matrix::colSums(x[["RNA"]]@counts[mito_genes, ])/total_counts_per_cell * 100

  # # Visualize QC metrics
  # VlnPlot(x, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  # Filter cells based on QC metrics
  # The thresholds here are just examples and might need adjustment based on your specific data
  # x <- subset(x, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  x <- subset(x, subset = nFeature_RNA > 200)

  return(x)
})


# ## Test seuratObject
# total_counts_per_cell <- colSums(seuratObject@assays$RNA@counts)
# mito_genes <- rownames(seuratObject)[grep("^MT-", rownames(seuratObject))]
# seuratObject$percent_mito <- colSums(seuratObject@assays$RNA@counts[mito_genes, ])/total_counts_per_cell
# VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
