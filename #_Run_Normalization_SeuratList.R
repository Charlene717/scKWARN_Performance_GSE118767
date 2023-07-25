##_Run_Normalization_SeuratList.R

## Set Normalization method
# Set_NorType <- c("counts","scKWARN", "RC", "scran", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
#                   "SCNorm", "sctransform", "PsiNorm")

## Call function
source("#_FUN_NorMeth.R")
seurat_Nor_list <- list()
logPlusOne <- Set_logPlusOne
# logPlusOne <- FALSE  # Set this variable to determine whether to apply log(x+1) transformation


# Apply function to decide whether to take log(x+1) of the data
applyLogPlusOne <- function(data, condition) {
  if(condition) {
    return(log(data + 1))
  } else {
    return(data)
  }
}

## Add to your source code
# Normalization methods as list of functions
normalization_methods <- list(
  "RC" = function(seuratObject, logPlusOne) {
    result <- naiveN(seuratObject@assays$RNA@counts)
    temp_object <- CreateSeuratObject(counts =seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(result$NormalizedData, logPlusOne) %>%  as(., "dgCMatrix")
    return(temp_object)
  },

  "Seurat_LogNormalize" = function(seuratObject, logPlusOne) {
    temp_object <- NormalizeData(seuratObject, normalization.method = "LogNormalize")
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(temp_object@assays[["RNA"]]@data, condition = FALSE)
    return(temp_object)
  },

  "Seurat_RC" = function(seuratObject, logPlusOne) {
    temp_object <- NormalizeData(seuratObject, normalization.method = "RC")
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(temp_object@assays[["RNA"]]@data, logPlusOne) %>%  as(., "dgCMatrix")
    return(temp_object)
  },

  "Seurat_CLR" = function(seuratObject, logPlusOne) {
    temp_object <- NormalizeData(seuratObject, normalization.method = "CLR")
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(temp_object@assays[["RNA"]]@data, logPlusOne)
    return(temp_object)
  },

  "scran" = function(seuratObject, logPlusOne) {
    result <- scranN2(seuratObject@assays$RNA@counts)
    temp_object <- CreateSeuratObject(counts = seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(result$NormalizedData, logPlusOne)%>%  as(., "dgCMatrix")
    return(temp_object)
  },

  "PsiNorm" = function(seuratObject, logPlusOne) {
    result <- PsiNorm(seuratObject  @assays$RNA@counts %>% as.matrix())
    temp_object <- CreateSeuratObject(counts = seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(result, logPlusOne) %>%  as(., "dgCMatrix")
    return(temp_object)
  },

  "scKWARN" = function(seuratObject, logPlusOne) {
    result <- LocASN(seuratObject@assays$RNA@counts, bw.method = "SJ",numGeneforEst = Set_numGeneforEst) # numGeneforEst = 500
    temp_object <- CreateSeuratObject(counts = seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(result$NormalizedData, logPlusOne) %>%  as(., "dgCMatrix")
    return(temp_object)
  },

  "SCNorm" = function(seuratObject, logPlusOne) {
    result <- scnormN(seuratObject@assays$RNA@counts, K = 5) #, K = 5 # Fix K if spent too much time
    temp_object <- CreateSeuratObject(counts =seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(result$NormalizedData, logPlusOne) %>%  as(., "dgCMatrix")
    return(temp_object)
  },


  "counts" = function(seuratObject, logPlusOne) {
    result <- seuratObject@assays$RNA@counts
    temp_object <- CreateSeuratObject(counts = seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    temp_object@assays[["RNA"]]@data <- applyLogPlusOne(result, logPlusOne) %>%  as(., "dgCMatrix")  # applyLogPlusOne(result, condition = FALSE)
    return(temp_object)
  },

  ## Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html
  "sctransform" = function(seuratObject, logPlusOne) {
    temp_object <- SCTransform(seuratObject, method = "nb")
    return(temp_object)
  },

  "sctransform_vst" = function(seuratObject, logPlusOne) {
    # Set_sctransform_Res = TRUE
    result <- sctrnN(seuratObject@assays$RNA@counts, setRes = Set_sctransform_Res)

    # Get the missing genes
    missing_genes <- setdiff(rownames(seuratObject@assays$RNA@counts), rownames(result$NormalizedData))

    # Fill missing values with 0
    missing_values <- matrix(0, nrow = length(missing_genes), ncol = ncol(result$NormalizedData))
    colnames(missing_values) <- colnames(result$NormalizedData)
    rownames(missing_values) <- missing_genes

    # Combine the normalized data and the missing values
    normalized_data <- rbind(result$NormalizedData, missing_values)

    # Reorder the rows to match the original SCE object
    normalized_data <- normalized_data[rownames(seuratObject@assays$RNA@counts), ]

    temp_object <- CreateSeuratObject(counts = seuratObject@assays[["RNA"]]@counts, meta.data = seuratObject@meta.data)
    if(Set_sctransform_Res == TRUE){
      temp_object@assays[["RNA"]]@data <- applyLogPlusOne(normalized_data, condition = FALSE) %>%  as(., "dgCMatrix")
    }else{
      temp_object@assays[["RNA"]]@data <- applyLogPlusOne(normalized_data, logPlusOne) %>%  as(., "dgCMatrix")
    }
    return(temp_object)
  }

  # Add more methods here...
)

## Apply normalization methods to each Seurat object
for (i in seq_along(seurat_list)) {
  # Get the current Seurat object
  seuratObject <- seurat_list[[i]]

  # Define a list to hold the normalized Seurat objects
  normalized_seurat_list <- list()

  # Apply each selected normalization method
  for (method in Set_NorType) {
    if (method %in% names(normalization_methods)) {
      try({
        normalized_seurat_list[[method]] <- normalization_methods[[method]](seuratObject, logPlusOne)
      })
    } else {
      print(paste("Normalization method", method, "not found. Skipping..."))
    }
  }

  # Save the list of normalized Seurat objects back into the original list
  seurat_Nor_list[[names(seurat_list)[i]]] <- normalized_seurat_list
}

rm(temp_object, seuratObject, result,
   missing_genes,missing_values)
