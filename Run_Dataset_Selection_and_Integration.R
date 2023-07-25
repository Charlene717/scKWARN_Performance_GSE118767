## Dataset Selection and Integration

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


