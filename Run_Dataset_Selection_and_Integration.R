## Ref: https://satijalab.org/seurat/articles/integration_introduction.html
## Ref: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)


##### Load required libraries #####
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)


# ##### Load data* #####
# seurat_list <- list()
# source("Dataset_PsiNorm.R")
#
# # source("Dataset_Seuratpbmc3k.R")
# # # source("Sup_Seurat_Sampling.R") # For Small test
# # seurat_list[["pbmc3k"]] <- seuratObject


#### Select Dataset ####
# Set_Dataset <- c( "mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53", "mix.DropSeq", "mix.10x",
#                   "mix.10x5", "mix.CELSeq" , "pbmc3k")

if(Name_DataSet=="mixCELSeq5"){
  # Define the datasets to include in the new seurat_list
  Set_Dataset <- c("mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53")

}else if(Name_DataSet=="mixCELSeq"){
  Set_Dataset <- c("mix.CELSeq")

}else if(Name_DataSet=="mix10x"){
  Set_Dataset <- c("mix.10x")

}else if(Name_DataSet=="mix10x5"){
  Set_Dataset <- c("mix.10x5")

}else if(Name_DataSet=="mixDropSeq"){
  Set_Dataset <- c("mix.DropSeq")

}else if(Name_DataSet=="pbmc3k"){
  Set_Dataset <- c("pbmc3k")

}else if(Name_DataSet==""){

}else{

}



# Remove datasets not in Set_Dataset from the seurat_list
seurat_list <- seurat_list[names(seurat_list) %in% Set_Dataset]


#### scRNA-seq integration ####
if(length(seurat_list) > 1){
  # Feature selection
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = Set_nfeatures)

  # Finding anchors
  seurat.anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)

  # Integrating data
  seurat.integrated <- IntegrateData(anchorset = seurat.anchors)

  # Storing the integrated data in the original objects
  seuratObject <- seurat.integrated

} else {
  seuratObject <- seurat_list[[1]]
}


# #### Perform an integrated analysis ####
# # Set_nfeatures = 2000
# if(length(seurat_list) > 1){
#   DefaultAssay(seuratObject) <- "integrated"
# } else {
#   seuratObject <- FindVariableFeatures(seuratObject, nfeatures = Set_nfeatures)
# }
#
#
# # Num_PCA = 50
# # Run the standard workflow for visualization and clustering
# seuratObject <- ScaleData(seuratObject, verbose = FALSE)
# seuratObject <- RunPCA(seuratObject, npcs = Num_PCA, verbose = FALSE)
# seuratObject <- RunUMAP(seuratObject, reduction = "pca", dims = 1:Num_PCA)
# seuratObject <- FindNeighbors(seuratObject, reduction = "pca", dims = 1:Num_PCA)
# seuratObject <- FindClusters(seuratObject, resolution = 0.5)
#
#
# # Visualization
# p1 <- DimPlot(seuratObject, reduction = "umap", group.by = "Dataset")
# p2 <- DimPlot(seuratObject, reduction = "umap", group.by = "seurat_clusters")
# # p2 <- DimPlot(seuratObject, reduction = "umap", label = TRUE, repel = TRUE)
# p3 <- DimPlot(seuratObject, reduction = "umap", group.by = "cell_line_demuxlet")
# p1 + p2 + p3
#
# ##### Count cell number #####
# source("#_FUN_Metric.R")
# Num_Cell.df <- calculate_counts(seuratObject, field_names = "cell_line_demuxlet")
# Num_Cell.df
#
