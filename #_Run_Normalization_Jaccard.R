## Perform normalization on the Jaccard program's data

## Call function
source("#_FUN_NorMeth.R")
Count.mtx <- seuratObject@assays[["RNA"]]@counts %>% as.matrix()

#### Normalizing the data ####
## RC ## scKWARN ## sctransform ## scran ## SCNorm ## PsiNorm
## scKWARN_SJ ## scKWARN_RoT ## sctransform ## scran ## Seurat_RC ## Seurat_CLR ## Seurat_LogNormalize

if(Name_NorType == "Seurat_LogNormalize"){
  ## Seurat_LogNormalize
  seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize")
  # seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

}else if(Name_NorType == "Seurat_CLR"){
  ## Seurat_CLR
  seuratObject <- NormalizeData(seuratObject, normalization.method = "CLR")

}else if(Name_NorType == "Seurat_RC"){
  ## Seurat_RC
  seuratObject <- NormalizeData(seuratObject, normalization.method = "RC")

}else if(Name_NorType == "RC"){
  Nor.mtx <- naiveN(Count.mtx)
  Nor.mtx <- as(Nor.mtx[["NormalizedData"]], "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

}else if(Name_NorType == "scKWARN"){
  ## scKWARN_SJ
  Nor.mtx <- LocASN(Count.mtx,bw.method = "SJ") #"RoT","SJ"
  Nor.mtx <- as(Nor.mtx[["NormalizedData"]], "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

}else if(Name_NorType == "SCNorm"){
  Nor.mtx <- scnormN(Count.mtx, K = 5) #, K = 5 # Fix K if spent too much time
  Nor.mtx <- as(Nor.mtx[["NormalizedData"]], "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

}else if(Name_NorType == "PsiNorm"){
  Nor.mtx <- PsiNorm(Count.mtx)
  Nor.mtx <- as(Nor.mtx, "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

}else if(Name_NorType == "scran"){
  ## scranN
  Nor.mtx <- scranN2(Count.mtx)
  Nor.mtx <- as(Nor.mtx[["NormalizedData"]], "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

## Ref: https://satijalab.org/seurat/articles/sctransform_vignette.html
}else if(Name_NorType == "sctransform"){
  ## sctransform
  seuratObject <- SCTransform(seuratObject, method = "nb")

}else if(Name_NorType == "sctransform_vst"){
  ## sctransform
  Nor.mtx <- sctrnN(Count.mtx, setRes = Set_sctransform_Res)
  seuratObject@assays[["RNA"]]@data <- Nor.mtx[["NormalizedData"]]

}else if(Name_NorType == "scKWARN_SJ"){
  ## scKWARN_SJ
  Nor.mtx <- LocASN(Count.mtx,bw.method = "SJ") #"RoT","SJ"
  # Nor.mtx <- scKWARN::LocASN(Count.mtx,bw.method = "SJ") #"RoT","SJ"
  Nor.mtx <- as(Nor.mtx[["NormalizedData"]], "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

}else if(Name_NorType == "scKWARN_RoT"){
  ## scKWARN_RoT
  Nor.mtx <- LocASN(Count.mtx,bw.method = "RoT") #"RoT","SJ"
  Nor.mtx <- as(Nor.mtx[["NormalizedData"]], "dgCMatrix")
  seuratObject@assays[["RNA"]]@data <- Nor.mtx

}

rm(Nor.mtx)

# # Test
# Nor.mtx_scranN2 <- seuratObject@assays[["RNA"]]@data
# Nor.mtx_Log <- seuratObject@assays[["RNA"]]@data
# sum(Nor.mtx_Log == Nor.mtx_scranN2)
