##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if (!requireNamespace("SeuratData", quietly = TRUE)) {install.packages("SeuratData")}; library(SeuratData)

cluster_Seurat_to_target <- function(seurat_obj, target_clusters = 5,
                                     Initial_resolution = 0.8, # Num_PCA = 50,
                                     LoopLim = 100, resolution_down_factor = 0.9,
                                     resolution_up_factor = 1.1) {

  num_clusters = 0
  counter = 1
  resolution <- Initial_resolution

  # Set a loop that runs at most LoopLim times to prevent infinite looping
  while (num_clusters != target_clusters && counter <= LoopLim) {
    # Perform clustering
    # seurat_obj <- FindNeighbors(seurat_obj, dims = 1:Num_PCA)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

    # Calculate the number of clusters
    num_clusters <- length(unique(Idents(seurat_obj)))

    # If there are too many clusters, lower the resolution parameter
    if (num_clusters > target_clusters) {
      resolution <- resolution * resolution_down_factor
    }

    # If there are too few clusters, raise the resolution parameter
    if (num_clusters < target_clusters) {
      resolution <- resolution * resolution_up_factor
    }

    # Update the counter
    counter <- counter + 1
  }

  # Store num_clusters and resolution to the seurat_obj
  seurat_obj@misc[["Rec_clustersNum"]]$num_clusters <- num_clusters
  seurat_obj@misc[["Rec_clustersNum"]]$resolution <- resolution
  seurat_obj@misc[["Rec_clustersNum"]]$target_clustersNum <- target_clusters

  # Return the seurat_obj
  return(seurat_obj)
}

# ## Test Function
# seurat_obj <- cluster_Seurat_to_target(seurat_obj, target_clusters = 5,
#                                        Initial_resolution = 0.5, Num_PCA = 50,
#                                        LoopLim = 100, resolution_down_factor = 0.9,
#                                        resolution_up_factor = 1.1)
#
# # Access the stored num_clusters and resolution
# print(paste("Number of clusters:", seurat_obj@misc[["Rec_clustersNum"]]$num_clusters))
# print(paste("Final resolution:", seurat_obj@misc[["Rec_clustersNum"]]$resolution))
#
#
# ## Run UMAP
# Num_PCA <- 50
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:Num_PCA)
# DimPlot(seurat_obj, reduction = "umap")
