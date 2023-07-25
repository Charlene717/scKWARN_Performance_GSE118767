## Dataset: GSE118767
## Ref: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767

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


# Assume that Dataset.list is a list of SingleCellExperiment objects
Dataset.list <- list(
  "mix.CELSeq51" = sc_Celseq2_5cl_p1,
  "mix.CELSeq52" = sc_Celseq2_5cl_p2,
  "mix.CELSeq53" = sc_Celseq2_5cl_p3,
  "mix.DropSeq" = sce_sc_Dropseq_qc,
  "mix.10x" = sce_sc_10x_qc,
  "mix.10x5" = sce_sc_10x_5cl_qc,
  # "pbmc3k" = pbmc3k,
  "mix.CELSeq" = sce_sc_CELseq2_qc
)

##### Integration #####
# Load required packages
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

# # Use the lapply function to convert each SingleCellExperiment object to Seurat object
# seurat_list <- lapply(Dataset.list, convert_to_seurat)
# Use the mapply function to convert each SingleCellExperiment object to Seurat object
seurat_list <- mapply(convert_to_seurat, Dataset.list, names(Dataset.list), SIMPLIFY=FALSE)


## Clean Up the Object
# rm(list = setdiff(ls(), c("Dataset.list", "seurat_list")))

# Get a list of all objects in the environment
objects <- ls()

# Find the objects that contain "sc_", "sce_", or "_qc"
objects_to_remove <- grep("sc_|sce_|_qc", objects, value = TRUE)

# Remove the selected objects
rm(list = objects_to_remove)
rm(objects,objects_to_remove)

