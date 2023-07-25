##### To-Do List ######
# - [V] Function of box plot


##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

## Record time set
Rec_Time_Point.lt <- list()
Rec_Time_Spend.lt <- list()

Rec_Time_Point.lt[["Start_Time"]] <- Sys.time() # %>% as.character()


##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)

if(!require("devtools")) install.packages("devtools")
devtools::install_github('satijalab/seurat-data'); library(SeuratData)
# if (!requireNamespace("SeuratData", quietly = TRUE)) {install.packages("SeuratData")}
# library(SeuratData)

## Install scKWARN by.tar.gz
library("scKWARN")

## Call function
source("#_FUN_NorMeth.R")
source("#_FUN_TimeRec.R")

Rec_Time_Point.lt[["Load_Packages"]] <- Sys.time() # %>% as.character()

##### Set Export #####
ExportFolder <- "Export"
if (!dir.exists(ExportFolder)){dir.create(ExportFolder)}   ## Create new folder

CMPName <- "MSINB"
DatasetName <- "SmallSimData" # "SmallSimData" # "Seurat_pbmc3k" # "GSE29087"
ExportName <- paste0(CMPName,"_",DatasetName,"_",
                     gsub("-", "_",Sys.Date()),"_",
                     "V1") # as.numeric(Sys.time(), units = "secs") %>% trunc())

##### Load data* #####
if(DatasetName == "Seurat_pbmc3k"){
  source("Dataset_Seuratpbmc3k.R")
  Count.mtx <- seuratObject@assays[["RNA"]]@counts %>% as.matrix()
  # Count.mtx <- as(Count.mtx, "sparseMatrix")

}else if(DatasetName == "GSE29087"){
  ## Load GSE29087 data
  source("Dataset_GSE29087_CountMtx.R")

}else if(DatasetName == "mix10x5"){
  source("Dataset_PsiNorm.R")
  Name_DataSet="mix10x5"
  source("Run_Dataset_Selection_and_Integration.R")
  seuratObject <- seurat_list[[1]]
  Count.mtx <- seuratObject@assays[["RNA"]]@counts %>% as.matrix()

}else{
  ## Simulation Example
  set.seed(12345)
  G <- 2000; n <- 600 # G: number of genes, n: number of cells
  NB_cell <- function(j) rnbinom(G, size = 0.1, mu = rgamma(G, shape = 2, rate = 2))
  countsimdata <- sapply(1:n, NB_cell)
  colnames(countsimdata) <- paste("cell", 1:n, sep = "_")
  rownames(countsimdata) <- paste("gene", 1:G, sep = "_")
  Count.mtx = countsimdata

}


Rec_Time_Point.lt[["Load_data"]] <- Sys.time() # %>% as.character()

##### Normalization #####
#### Seurat4 ####
# Ref: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

Results_Time.df <- data.frame(Method = character(),Times = integer(), TimeElapsed = numeric())

methods_list <- list(
  list(Method_Name = "scKWARN", func = scKWARN::LocASN, bw.method = "SJ"),
  # list(Method_Name = "scKWARN_RoT", func = scKWARN::LocASN, bw.method = "RoT"),
  # list(Method_Name = "scKWARN_SJ", func = scKWARN::LocASN, bw.method = "SJ"),
  # list(Method_Name = "Seurat_Log", func = Seurat::NormalizeData, normalization.method = "LogNormalize"),
  # list(Method_Name = "Seurat_CLR", func = Seurat::NormalizeData, normalization.method = "CLR"),
  # list(Method_Name = "Seurat_RC", func = Seurat::NormalizeData, normalization.method = "RC"),
  list(Method_Name = "RC", func = naiveN),
  list(Method_Name = "scran", func = scranN2),
  list(Method_Name = "SCNorm", func = scnormN, K = 5),
  list(Method_Name = "sctransform", func = sctrnN),
  list(Method_Name = "PsiNorm", func = PsiNorm)
)

default_args <- list(InputData = Count.mtx, Set_Times = 10)

for (method in methods_list) {
  Test.lt <- do.call(FUN_TimeRec, modifyList(default_args, modifyList(method, list(TimePoint.lt = Rec_Time_Point.lt, TimeSpend.lt = Rec_Time_Spend.lt, Time.df = Results_Time.df))))
  Rec_Time_Point.lt <- Test.lt[["TimePoint.lt"]]
  Rec_Time_Spend.lt <- Test.lt[["TimeSpend.lt"]]
  Results_Time.df <- Test.lt[["Time.df"]]
  rm(Test.lt)
}


##### Record time df #####
## Record time df
Results_Time.df <- data.frame(ID = ExportName,
                              CMPName = CMPName,
                              DatasetName = DatasetName,
                              CellNum = ncol(Count.mtx),
                              GeneNum = nrow(Count.mtx),
                              Results_Time.df
)


##### Plot #####
## ggplot2 box plot
## Ref: http://www.sthda.com/english/wiki/ggplot2-box-plot-quick-start-guide-r-software-and-data-visualization

# ## Load data
# ImportName <- "MSINB_0511_0513_Test/MSINB_SmallSimData_2023_05_11_V1"
# Results_Time.df <- read.delim2(paste0(getwd(),"/Export/",ImportName,"_TimeRec.tsv"))
# load(paste0(getwd(),"/Export/",ImportName,".RData"))

source("FUN_Plot_Box.R")
# Set Color
col.values <- c("counts" = "gray","scKWARN" = "blue", "RC" = "black", "scran" = "red", "SCNorm" = "#00CC00",
                "Seurat_LogNormalize" = "#159e5a","Seurat_RC" = "#2bed78","Seurat_CLR" = "#335e49",
                "sctransform" = "#B266FF", "PsiNorm" = "#FF8000")
Set_NorType <- c("counts","scKWARN", "RC", "scran", "SCNorm",  # "SCNorm", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
                 "sctransform", "PsiNorm")

NumCell <- ncol(Count.mtx)
NumGene <- nrow(Count.mtx)

plt.box <- create_box_plot(Results_Time.df,
                           x_var = "Method", y_var = "TimeElapsed",
                           fill_var = "Method", color_var = "Method",
                           Set_order = "user",
                           Name_XTitle = "Method", Name_YTitle = "Time spent(Sec)",
                           Set_Title = paste0(DatasetName," (",NumCell,"cells, ",NumGene,"genes)"),
                           color_values = col.values,
                           level_values = Set_NorType)
plt.box


##### Export #####
write.table(Results_Time.df,
            file=paste0(ExportFolder,"/",ExportName,"_TimeRec.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')
## PDF
ggsave(filename = paste0(ExportFolder,"/",ExportName,"_BarPlot.pdf"), plot = plt.box, width = 8, height = 8, units = "in")
## TIFF
ggsave(filename = paste0(ExportFolder,"/",ExportName,"_BarPlot.tiff"), plot = plt.box, width = 8, height = 8, units = "in")

## Record version and sessionInfo
info_output <- c("##_R Version Information:", capture.output(version), "",
                 "##_Session Information:", capture.output(sessionInfo()))

writeLines(info_output, paste0(ExportFolder,"/",ExportName,"_Version_and_Session_Info.txt"))

## RData
save.image(paste0(ExportFolder,"/",ExportName,"_TimeRec.RData"))

