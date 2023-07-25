# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)

##### Load data* #####
# ## Seuratpbmc3k backup method (Load from folder)
# Count.mtx <- read.delim2("./Input_Count/SeuratPBMC_CountMTX.tsv") %>%
#   magrittr::set_rownames(.[, 1]) %>%
#   dplyr::select(-1) %>%
#   as.matrix()

## Seuratpbmc3k from SeuratData
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

# AvailableData()
InstallData("pbmc3k")
data("pbmc3k")
seuratObject <- pbmc3k
rm(pbmc3k)

# Define a function to deal with missing values in the metadata of a Seurat object
handle_na_seurat <- function(seurat_object, column_name, replacement_value = NA, remove_na = FALSE) {
  # Find rows with missing values or empty strings in the metadata
  is_missing <- is.na(seurat_object@meta.data[, column_name]) | seurat_object@meta.data[, column_name] == ""
  missing_rows <- which(is_missing)

  if (length(missing_rows) > 0) {
    cat("Rows with missing values:", missing_rows, "\n")

    if (remove_na) {
      # Remove rows with missing values
      seurat_object@meta.data <- seurat_object@meta.data[-missing_rows, ]
      seurat_object@assays$RNA@data <- seurat_object@assays$RNA@data[, -missing_rows]
    } else {
      # Replace missing values with the specified replacement value
      column_data <- as.character(seurat_object@meta.data[, column_name])
      column_data[is_missing] <- replacement_value
      seurat_object@meta.data[, column_name] <- as.factor(column_data)
    }
  }

  return(seurat_object)
}



# Use the function to replace missing values with "Unknown"
seuratObject <- handle_na_seurat(seuratObject, 'seurat_annotations', "Unknown")

# # Use the function to remove rows with missing values
# seuratObject <- handle_na_seurat(seuratObject, 'seurat_annotations', remove_na = TRUE)


# ##### Check results #####
# # Run the standard workflow for visualization and clustering
# seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
# seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
# seuratObject <- ScaleData(seuratObject, verbose = FALSE)
# seuratObject <- RunPCA(seuratObject, npcs = 30, verbose = FALSE)
# seuratObject <- RunUMAP(seuratObject, reduction = "pca", dims = 1:30)
# seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:30)
# seuratObject <- FindClusters(seuratObject, resolution = 0.5)
# # Visualization
# p1 <- DimPlot(seuratObject, reduction = "umap", group.by = "seurat_annotations")
# p2 <- DimPlot(seuratObject, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2
#
# rm(p1,p2)
