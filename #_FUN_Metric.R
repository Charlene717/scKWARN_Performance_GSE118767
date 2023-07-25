##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("cluster")) install.packages("cluster"); library(cluster)
if(!require("mclust")) install.packages("mclust"); library(mclust)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if(!require("ROGUE")) devtools::install_github("PaulingLiu/ROGUE")


##### Metric #####
#### PCADepthCorr ####
calculate_PCADepthCorr <- function(seurat_object, numPCs = 2, logDepth = FALSE) {
  # Initialize the PCADepthCorr list
  PCADepthCorr <- list()
  PCADepthCorr$max_correlation <- NA
  max_corr_PC <- NA

  # Calculate the depth
  depth <- colSums(seurat_object@assays[["RNA"]]@counts)

  # Apply Log transformation on depth if logDepth is TRUE
  if (logDepth) {
    depth <- log(depth)
  }

  # Calculate the correlation with sequencing depth for each PC
  for (i in 1:numPCs) {
    PC <- seurat_object@reductions$pca@cell.embeddings[, i]
    corr_PC <- abs(cor(PC, depth))

    # Add the correlation to the PCADepthCorr list
    PCADepthCorr[[paste0("PC", i)]] <- round(corr_PC, digits = 3)

    # Update the maximum correlation
    if (is.na(PCADepthCorr$max_correlation) || corr_PC > PCADepthCorr$max_correlation) {
      PCADepthCorr$max_correlation <- round(corr_PC, digits = 3)
      max_corr_PC <- paste0("PC", i)
    }
  }

  PCADepthCorr$max_correlation_PC <- max_corr_PC

  return(PCADepthCorr)
}



#### Average Silhouette Width ####
# Load the necessary libraries
library(cluster)
library(Seurat)

calculate_silhouette <- function(seurat_object, reduction_method = "umap",
                                 distance_method = "euclidean", cluster_method = "seurat_clusters",
                                 pca_components = NULL) {

  # Check if reduction_method is "pca" and pca_components is provided
  if (reduction_method == "pca" && is.null(pca_components)) {
    stop("Number of PCA components must be specified when reduction_method is 'pca'")
  }

  # Extract the cell embeddings based on the reduction_method
  if (reduction_method == "pca") {
    cell_embeddings <- Embeddings(object = seurat_object[["pca"]])[, 1:pca_components]

  } else {
    cell_embeddings <- Embeddings(seurat_object, reduction = reduction_method)
  }

  # Extract the cluster assignments
  cluster_assignments <- seurat_object@meta.data[[cluster_method]]

  # If cluster assignments are factor or character type, convert them to numeric labels
  if (is.factor(cluster_assignments) || is.character(cluster_assignments)) {
    cluster_assignments <- as.integer(as.factor(cluster_assignments))
  }

  # Ensure cluster assignments are numeric and not NA
  if (any(is.na(cluster_assignments))) stop("NA values found in 'cluster_assignments'")

  # Check if the cluster_assignments vector contains only integers
  if (!all(cluster_assignments == round(cluster_assignments))) stop("'cluster_assignments' must contain only integer values")

  # Calculate the dissimilarity matrix using the specified distance method
  dissimilarity_matrix <- dist(cell_embeddings, method = distance_method)

  # Calculate the silhouette values
  silhouette_values <- silhouette(cluster_assignments, dissimilarity_matrix)

  # Calculate the average silhouette width
  average_silhouette_width <- mean(silhouette_values[, 3])

  # Return a list that includes the updated object and the silhouette values
  return(list("silhouette_values" = silhouette_values,
              "average_silhouette_width" = average_silhouette_width))
}


#### Purity ####
calculate_purity <- function(seurat_object,
                             annotation_name = "seurat_annotations",
                             cluster_name = "seurat_clusters") {
  # Extract cell types and clusters
  cell_types <- seurat_object@meta.data[[annotation_name]]
  clusters <- seurat_object@meta.data[[cluster_name]]

  # Compute the distribution of cell types in each cluster
  cell_type_distribution <- table(cell_types, clusters)

  # Compute the purity of each cluster (proportion of the most common cell type)
  cluster_purity <- apply(cell_type_distribution, 2, function(x) max(x) / sum(x))

  # Compute the most common cell type in each cluster
  main_cell_types <- rownames(cell_type_distribution)[apply(cell_type_distribution, 2, function(x) which.max(x))]

  # Create a data frame that includes the cluster names, purity, and most common cell type
  cluster_purity_df <- data.frame(cluster_id = names(cluster_purity),
                                  Purity = as.numeric(cluster_purity),
                                  Main_Cell_Type = main_cell_types,
                                  stringsAsFactors = FALSE)

  # Return the data frame
  return(cluster_purity_df)
}


#### ROGUE ####
## Ref: https://htmlpreview.github.io/?https://github.com/PaulingLiu/ROGUE/blob/master/vignettes/ROGUE_Tutorials.html

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if(!require("ROGUE")) devtools::install_github("PaulingLiu/ROGUE")

# Load the necessary library
library(tidyverse)
library(dplyr)
library(ROGUE) # suppressMessages(library(ROGUE))

calculate_rogue_value <- function(seuratObject, labels, samples, min.cells = 10, min.genes = 10){
  # Load the data
  expr <- seuratObject@assays[["RNA"]]@data %>% as.matrix()
  meta <- seuratObject@meta.data

  # Ensure labels and samples are in meta data
  if(!(labels %in% names(meta))){
    stop(paste0("The label '", labels, "' is not found in the meta data of the Seurat object."))
  }

  if(!(samples %in% names(meta))){
    stop(paste0("The sample '", samples, "' is not found in the meta data of the Seurat object."))
  }

  ## Filtering out low-abundance genes and low-quality cells
  expr <- matr.filter(expr, min.cells = min.cells, min.genes = min.genes)

  ## Expression entropy model
  ent.res <- SE_fun(expr)

  ## Calculate the ROGUE value of each putative cluster for each sample
  rogue.res <- rogue(expr, labels = meta[[labels]], samples = meta[[samples]], platform = "UMI", span = 0.6)

  # Visualize ROGUE values on a boxplot
  rogue.boxplot(rogue.res)

  # Calculate mean rogue value, excluding NAs
  avg_rogue.value <- mean(unlist(rogue.res), na.rm = TRUE)

  ## Output
  Output <- list(avg_rogue.value = avg_rogue.value,
                 rogue.df = rogue.res)

  return(Output)
}



#### Adjusted Rand index (ARI) ####
# Load the necessary libraries
library(Seurat)
library(mclust)

calculate_ari <- function(seurat_obj, true_labels_field, predicted_labels_field) {
  # Get the true cell types
  true_labels <- as.factor(seurat_obj@meta.data[[true_labels_field]])

  # Get the predicted cell types
  predicted_labels <- as.factor(seurat_obj@meta.data[[predicted_labels_field]])

  # Calculate the ARI using mclust's adjustedRandIndex function
  ari <- adjustedRandIndex(true_labels, predicted_labels)

  # Return the ARI
  return(ari)
}

# ## Test function
# ari <- calculate_ari(seurat_obj, "Cell_Type", "seurat_clusters")
# print(ari)










####################################################################################################
#### Count cell number ####
calculate_counts <- function(seurat_object, field_names) {
  # Initialize a list to store the data frames
  df_list <- list()

  # Loop over each field name
  for (field_name in field_names) {
    # Extract the specified field from the Seurat object
    field_data <- seurat_object@meta.data[[field_name]]

    # Compute the count of each category in the field
    field_counts <- table(field_data)

    # Convert the categories and their counts to a data frame
    # Change the column name 'ID' to the value of field_name
    df <- data.frame(field_name = names(field_counts),
                     count = as.vector(field_counts))
    names(df)[1] <- field_name

    # Add the data frame to the list
    df_list[[field_name]] <- df
  }

  # Return the list of data frames
  return(df_list)
}



