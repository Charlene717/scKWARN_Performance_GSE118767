#### Run A549 ####
rm(list = ls())
memory.limit(150000)

Set_Test_Type <- "CompDiff"
Set_Main_CellLine <- "A549"
Set_Main_FltCellLine = "ALL"
Set_TotalCellNum <- 800
Set_numGeneforEst = 2000
Set_Filter_Cell = TRUE

Set_NorType <- c("counts","scKWARN", "RC", "scran", "sctransform", "PsiNorm")
to_keep <- c("i", "Name_Test", "Set_Test_Type", "Set_Main_CellLine", "Set_TotalCellNum" ,
             "Set_NorType","Set_numGeneforEst", "Set_Filter_Cell","Set_Main_FltCellLine",
             "Set_FindVarFea", "varFea_values", "to_keep")

# Here we create a vector with the different values for Set_FindVarFea
varFea_values <- "vst" # varFea_values <- c("disp", "mvp", "vst")

for(Set_FindVarFea in varFea_values) {
  for(i in 1:20) {
    # set.seed(123+i)
    Set_Seed <- i+123
    Name_Test <- paste0("AllF_V", i)

    source("##_scKWARN_Test_RealData_IntegrAll_ForRunAll_20230625.R")

    # Remove all objects except those in to_keep
    rm(list=setdiff(ls(), to_keep))

    # Collect unused memory
    gc()
  }
}



#### Run H838 ####
rm(list = ls())
memory.limit(150000)

Set_Test_Type <- "CompDiff"
Set_Main_CellLine <- "H838"
Set_Main_FltCellLine = "ALL"
Set_TotalCellNum <- 500
Set_numGeneforEst = 2000
Set_Filter_Cell = TRUE

Set_NorType <- c("counts","scKWARN", "RC", "scran", "sctransform", "PsiNorm")
to_keep <- c("i", "Name_Test", "Set_Test_Type", "Set_Main_CellLine", "Set_TotalCellNum" ,
             "Set_NorType","Set_numGeneforEst", "Set_Filter_Cell","Set_Main_FltCellLine",
             "Set_FindVarFea", "varFea_values", "to_keep")

# Here we create a vector with the different values for Set_FindVarFea
varFea_values <- "vst" # varFea_values <- c("disp", "mvp", "vst")

for(Set_FindVarFea in varFea_values) {
  for(i in 1:20) {
    # set.seed(123+i)
    Set_Seed <- i+123
    Name_Test <- paste0("AllF_V", i)

    source("##_scKWARN_Test_RealData_IntegrAll_ForRunAll_20230625.R")

    # Remove all objects except those in to_keep
    rm(list=setdiff(ls(), to_keep))

    # Collect unused memory
    gc()
  }
}

