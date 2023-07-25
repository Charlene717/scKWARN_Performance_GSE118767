## Function for Normalization

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

## scran
if(!require("scran")) BiocManager::install("scran", dependencies = T); library(scran)

scranN2 <- function(counts2, clusters = NULL) {
  sce <- SingleCellExperiment(assays=list(counts=counts2))
  if(!is.null(clusters)) {
    clusters <- quickCluster(counts2)
    sce <- computeSumFactors(sce, clusters = clusters)$sizeFactor
  } else {
    sce <- computeSumFactors(sce, clusters = clusters)$sizeFactor
  }
  list(NormalizedData =  t(t(counts2)/sce), scalingFactor = sce)
}


## SCNorm
if(!require("SCnorm")) BiocManager::install("SCnorm", dependencies = T); library(SCnorm)
scnormN <- function(counts, conditions = rep(1, ncol(counts)), K = NULL) {
  DataNorm <- SCnorm(Data = data.frame(counts),
                     Conditions = conditions,
                     PrintProgressPlots = F,
                     FilterCellNum = 10, K = K,
                     NCores = 2, reportSF = TRUE)
  list(NormalizedData =  SingleCellExperiment::normcounts(DataNorm))
}


## sctransform
if(!require("sctransform")) BiocManager::install("sctransform", dependencies = T); library(sctransform)
sctrnN <- function(counts,setRes = FALSE) {
  sparse_data <- as(as.matrix(counts), "sparseMatrix")
  sct_results <- sctransform::vst(umi = sparse_data, return_cell_attr = TRUE, return_gene_attr = TRUE, return_corrected_umi = TRUE, method = "nb")
  if(setRes == TRUE){
    SCT_data <- sct_results$y # Residual
  }else{
    SCT_data <- sct_results$umi_corrected
  }
  list(NormalizedData = SCT_data)
}

# # Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html
# seuratObject <- SCTransform(seuratObject)

## RC LibrarySize
naiveN <- function(counts) {
  sj <- colSums(counts)
  mg   <- sj/median(sj)
  SN <- t(t(counts)/mg)
  SN <- as(SN, "sparseMatrix")
  list(NormalizedData = SN, scalingFactor = mg)
}
## Test function
# M <- matrix(c(1,0,2,0,2,9,3,0), ncol=2)
# naiveN(M)


## PsiNorm
# Ref:https://pubmed.ncbi.nlm.nih.gov/34499096/
if(!require("PsiNorm"))  remotes::install_github("MatteoBlla/PsiNorm"); library(PsiNorm)


## scKWARN
## Install scKWARN by.tar.gz
library("scKWARN")
# ## Backup scKWARN loading method
# source("LocASN.R")
# source("Core_functions.R")
# library(Rcpp)
# sourceCpp("calculateLocAve.cpp")
# Nor.mtx <- LocASN(Count.mtx)


