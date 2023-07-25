# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)

# Load required libraries
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require("SingleCellExperiment")) BiocManager::install("SingleCellExperiment"); library(SingleCellExperiment)
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)


##### Load Dataset #####
# Automatically get current working directory and append your subfolder
dir_path <- paste0(getwd(), "/Input_Dataset/GSE118767/")

# List all .RData files in the directory
files <- list.files(path = dir_path, pattern = "*.RData", full.names = TRUE)

# Load all .RData files
lapply(files, load, .GlobalEnv)


mix.10x5 <- sce_sc_10x_5cl_qc

#### convert SingleCellExperiment to Seurat ####
library(Seurat)
library(SingleCellExperiment)

# Define a function to convert SingleCellExperiment to Seurat
convert_to_seurat <- function(sce, dataset_name) {

  # Ensure that we are working with counts data rather than logcounts
  counts <- counts(sce)

  # Create a new Seurat object
  seurat_object <- CreateSeuratObject(counts = counts)

  # # Copy the metadata from the SingleCellExperiment object to the Seurat object
  # metadata(sce) <- seurat_object@meta.data

  # Copy the metadata from the SingleCellExperiment object to the Seurat object
  seurat_object@meta.data <- data.frame(seurat_object@meta.data, colData(sce))

  # Add a "Dataset" field in metadata, and set it to the dataset_name
  seurat_object$Dataset <- dataset_name


  return(seurat_object)
}

seurat_object <- convert_to_seurat(mix.10x5,"mix.10x5")

## Clean Up the Object
# Get a list of all objects in the environment
objects <- ls()

# Find the objects that contain "sc_", "sce_", or "_qc"
objects_to_remove <- grep("sc_|sce_|_qc", objects, value = TRUE)

# Remove the selected objects
rm(list = objects_to_remove)
rm(mix.10x5,objects,objects_to_remove)

############################################################################################
## Apply the preprocessing steps for Seurat Object
# Set_nfeatures <- 2000
# Num_PCA <- 50
seurat_object <- seurat_object %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = Set_nfeatures) %>%
  ScaleData() %>%
  RunPCA(npcs = Num_PCA) %>%
  FindNeighbors(dims = 1:Num_PCA) %>%
  FindClusters(resolution = 0.5) %>% # resolution = 0.5
  RunUMAP(dims = 1:Num_PCA) %>%
  RunTSNE(dims = 1:Num_PCA)

DimPlot(seurat_object, reduction = "umap",label = TRUE)
DimPlot(seurat_object, reduction = "umap", group.by = "cell_line_demuxlet",label.size = 2)
DimPlot(seurat_object, reduction = "umap", group.by = "cell_line",label.size = 2)
DimPlot(seurat_object, reduction = "tsne", group.by = "cell_line_demuxlet")

## Clean up the cluster
# Suppose the name of your Seurat object is seurat_object
# First, determine what your cluster labels are, usually stored in the "seurat_clusters" field of metadata
table(seurat_object$seurat_clusters) # List all clusters, from which you can see the label of the 5th cluster

# Use the subset function to retain all clusters except specific cluster
if (Set_Main_FltCellLine == "A549") {
  seurat_object_filtered <- seurat_object[,which(!seurat_object$seurat_clusters %in% c(5))]
} else if (Set_Main_FltCellLine == "H838") {
  seurat_object_filtered <- seurat_object[,which(!seurat_object$seurat_clusters %in% c(6,7))]
} else {
  seurat_object_filtered <- seurat_object[,which(!seurat_object$seurat_clusters %in% c(5,6,7,8,9,10))]
}

Plt_UMAP_Flt <- DimPlot(seurat_object_filtered, reduction = "umap",label = TRUE) +
  theme(panel.background = element_rect(fill = NA, color = "black", size = 1.2)) +
  coord_fixed()
Plt_UMAP_Ori <-DimPlot(seurat_object, reduction = "umap",label = TRUE) +
  theme(panel.background = element_rect(fill = NA, color = "black", size = 1.2)) +
  coord_fixed()
Plt_UMAP <- Plt_UMAP_Ori + Plt_UMAP_Flt
Plt_UMAP
Plt_UMAP_CT_Ori <- DimPlot(seurat_object, reduction = "umap", group.by = "cell_line_demuxlet",label.size = 2)+
  theme(panel.background = element_rect(fill = NA, color = "black", size = 1.2)) +
  coord_fixed()
Plt_UMAP_CT_Flt <- DimPlot(seurat_object_filtered, reduction = "umap", group.by = "cell_line_demuxlet",label.size = 2)+
  theme(panel.background = element_rect(fill = NA, color = "black", size = 1.2))+
  coord_fixed()
Plt_UMAP_CT <- Plt_UMAP_CT_Ori + Plt_UMAP_CT_Flt
Plt_UMAP_CT

seurat_F_list <- list()
seurat_F_list[["mix.10x5"]] <- seurat_object
seurat_F_list[["mix.10x5_filtered"]] <- seurat_object_filtered

## Plot CellCount
source("#_FUN_Metric.R")
seurat_F_list <- lapply(seurat_F_list, function(seurat_object) {
  # Calculate counts
  counts <- calculate_counts(seurat_object, field_names = "cell_line_demuxlet")

  # Store the result in the misc slot
  seurat_object@misc[["CountCell"]] <- counts

  return(seurat_object)
})

## Plot Result
source("FUN_create_plots_CellCount_list.R")
Set_YLimCellNum <- 1350

plots_CellCount_list <- create_plots_CellCount_list(seurat_F_list,Set_Cell = "cell_line_demuxlet",YLimCellNum =Set_YLimCellNum)
grid.arrange(grobs = plots_CellCount_list, ncol = 3)

## Export PDF
pdf(paste0(Name_ExportFolder,"/",Name_Export,"_FilterCell.pdf"),
    width = 16, height = 8)
grid.arrange(grobs = plots_CellCount_list, ncol = 3)
print(Plt_UMAP)
print(Plt_UMAP_CT)
dev.off()

#### Run ROGUE ####
## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html
## Load package
suppressMessages(library(ROGUE))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))

seurat_object <- seurat_F_list[["mix.10x5"]]
expr <- seurat_object@assays[["RNA"]]@data %>% as.matrix()
meta <- seurat_object@meta.data %>% as.data.frame()

## Filtering out low-abundance genes and low-quality cells
expr <- matr.filter(expr, min.cells = 10, min.genes = 10)
## Expression entropy model
ent.res <- SE_fun(expr)
head(ent.res)
SEplot(ent.res)

## ROGUE calculation
rogue.value <- CalculateRogue(ent.res, platform = "UMI")
rogue.value

rogue.res <- rogue(expr, labels = meta$cell_line_demuxlet, samples = meta$cell_line_demuxlet, platform = "UMI", span = 0.6)
rogue.res

#Visualize ROGUE values on a boxplot
rogue.boxplot(rogue.res)



