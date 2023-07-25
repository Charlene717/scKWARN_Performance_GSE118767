## Main Code: Main Metrics (Run All situation)


#### Run GSE118767 ####
rm(list = ls())
memory.limit(150000)

Set_Test_Type <- ""
# Set_Main_CellLine <- "A549"
# Set_Filter_Cell = FALSE
# Set_Main_FltCellLine = "ALL"
# Set_TotalCellNum <- 800
Set_numGeneforEst = 2000
Set_FindVarFea <- "vst"
Set_NorType <- c("counts","scKWARN", "RC", "scran", "sctransform", "PsiNorm", "SCNorm")

source("##_scKWARN_RealData_IntegrAll.R")

##################################################################################

##### RunLoop: Compositional Differences #####
#### Run A549 ####
rm(list = ls())
memory.limit(150000)

Set_Test_Type <- "CompDiff"
Set_Main_CellLine <- "A549"
Set_Filter_Cell = TRUE
Set_Main_FltCellLine = "ALL"
Set_TotalCellNum <- 800
Set_numGeneforEst = 2000
varFea_values <- "vst" # Here we create a vector with the different values for Set_FindVarFea

Set_NorType <- c("counts","scKWARN", "RC", "scran", "sctransform", "PsiNorm")
to_keep <- c("i", "Name_Test", "Set_Test_Type", "Set_Main_CellLine", "Set_TotalCellNum" ,
             "Set_NorType","Set_numGeneforEst", "Set_Filter_Cell","Set_Main_FltCellLine",
             "Set_FindVarFea", "varFea_values", "to_keep")

for(Set_FindVarFea in varFea_values) {
  for(i in 1:20) {
    # set.seed(123+i)
    Set_Seed <- i+123
    Name_Test <- paste0("AllF_V", i)

    source("##_scKWARN_RealData_IntegrAll.R")

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
Set_Filter_Cell = TRUE
Set_Main_FltCellLine = "ALL"
Set_TotalCellNum <- 500
Set_numGeneforEst = 2000
varFea_values <- "vst" # Here we create a vector with the different values for Set_FindVarFea

Set_NorType <- c("counts","scKWARN", "RC", "scran", "sctransform", "PsiNorm")
to_keep <- c("i", "Name_Test", "Set_Test_Type", "Set_Main_CellLine", "Set_TotalCellNum" ,
             "Set_NorType","Set_numGeneforEst", "Set_Filter_Cell","Set_Main_FltCellLine",
             "Set_FindVarFea", "varFea_values", "to_keep")

for(Set_FindVarFea in varFea_values) {
  for(i in 1:20) {
    # set.seed(123+i)
    Set_Seed <- i+123
    Name_Test <- paste0("AllF_V", i)

    source("##_scKWARN_RealData_IntegrAll.R")

    # Remove all objects except those in to_keep
    rm(list=setdiff(ls(), to_keep))

    # Collect unused memory
    gc()
  }
}


##################################################################################

#### Visualization: Box Plot ####
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

##### Set parameter* ####
Set_FindVarFea <- "vst" # "disp" # "mvp" # "vst"

Name_Title = paste0("CompDiff_mix10x5_A549_Sum800 nfeature:2000; numGeneforEst:2000; nPC:50; FindVarFea:",Set_FindVarFea)
# Name_Title = paste0("CompDiff_mix10x5_H838_Sum500 nfeature:2000; numGeneforEst:2000; nPC:50; FindVarFea:",Set_FindVarFea)

Name_ExportFolder = paste0("Export/#_MSINB_0625_Performance_CompDiff_A549_scKW2000_FilterAllCell_#")
# Name_ExportFolder = paste0("Export/#_MSINB_0625_Performance_CompDiff_H838_scKW2000_FilterAllCell_#")

Name_Export = paste0(gsub("-", "_",Sys.Date()),"_Sum_Metics_", Set_FindVarFea)
Set_SelectPlt <- c("ASWPCA_CellType", "ARI", "Inverse_PCADepthCorr" )

##### Load data #####
Set_Path <- paste0(getwd(), "/" ,Name_ExportFolder)

## Merge specific files from a folder into a Dataframe
files <- list.files(Set_Path, pattern = ".tsv", full.names = TRUE)
files <- list.files(Set_Path, pattern = "*Sum*.tsv", full.names = TRUE)

## Read each tsv file and stack them into a data frame
Sum.df <- files %>%
  map_dfr(read_tsv)

print(Sum.df) # Display the data

# Check if the column "PCADepthCorr" exists
if ("PCADepthCorr" %in% colnames(Sum.df)) {
  # If exists, create a new column "Inverse_PCADepthCorr"
  Sum.df$Inverse_PCADepthCorr <- 1 - Sum.df$PCADepthCorr
}

## Filter some condition
Sum.df <- Sum.df[!Sum.df$NorMeth %in% "counts" & !Sum.df$Dataset %in% "Sum_800_A549_80", ]


##### Plot Box #####
library(ggplot2)
source("FUN_Plot_Box.R")
source("FUN_Plot_Arrange.R")

# Set Color
col.values <- c("counts" = "gray","scKWARN" = "blue", "RC" = "black", "scran" = "red", "SCNorm" = "#00CC00",
                "Seurat_LogNormalize" = "#159e5a","Seurat_RC" = "#2bed78","Seurat_CLR" = "#335e49",
                "sctransform" = "#B266FF", "PsiNorm" = "#FF8000")
Set_NorType <- c("counts","scKWARN", "RC", "scran",  # "SCNorm", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
                 "sctransform", "PsiNorm")
# Set_Dataset <- c("mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53", "mix.DropSeq", "mix.10x", "mix.10x5", "mix.CELSeq")
Set_Dataset <- Sum.df[,"Dataset"] %>% unique() %>% unlist() %>% rev()

# Define fixed parameters for create_box_plot_multiXTpye function
fixed_params <- list(x_var = "Dataset", fill_var = "NorMeth",color_var = "NorMeth",
                     color_values = col.values,show_mean = FALSE, # ylim_values = c(0, 1),
                     Set_fill_order = Set_NorType, Set_x_order = Set_Dataset,
                     box_alpha = 0.1, jitter_alpha = 0.7, box_line = 0.8, AspectRatio = 0.8,
                     x_text_size = 17, y_text_size = 17, legend_text_size = 26,
                     x_title_size = 32, y_title_size = 32,legend_title_size = 24,
                     Name_XTitle = "Dataset", Name_LegendTitle = "",
                     legend_direction = "horizontal")
new_params <- list(AspectRatio = 1, show_mean = TRUE, ylim_values = c(0, 1), rect_size = 2)
fixed_params2 <- modifyList(fixed_params, new_params)


# Define list for variable parameters with names
variable_params_list <- list(
  PCADepthCorr = list(df = Sum.df, y_var = "PCADepthCorr", Name_YTitle = "Correlation"),
  Inverse_PCADepthCorr = list(df = Sum.df, y_var = "Inverse_PCADepthCorr", Name_YTitle = "1-Correlation"),
  ASWPCA_Cluster = list(df = Sum.df, y_var = "ASWPCA_Cluster", Name_YTitle = "ASW (Cluster)"),
  ASWPCA_CellType = list(df = Sum.df, y_var = "ASWPCA_CellType", Name_YTitle = "ASW"),
  ARI = list(df = Sum.df, y_var = "ARI", Name_YTitle = "ARI"),
  Purity = list(df = Sum.df, y_var = "Purity", Name_YTitle = "Purity"),
  ROGUE_Cluster = list(df = Sum.df, y_var = "ROGUE_Cluster", Name_YTitle = "ROGUE (Cluster)"),
  ROGUE_CellType = list(df = Sum.df, y_var = "ROGUE_CellType", Name_YTitle = "ROGUE (CellType)")
)

variable_params_list <- variable_params_list[Set_SelectPlt]

# Use lapply to apply create_box_plot_multiXTpye function
plot_list_Box <- lapply(variable_params_list, function(params) {
  args <- modifyList(fixed_params, params)
  do.call(create_box_plot_multiXTpye, args)
})
plot_list_Box2 <- lapply(variable_params_list, function(params) {
  args <- modifyList(fixed_params2, params)
  do.call(create_box_plot_multiXTpye, args)
})

Name_FigTitle <- Name_Title
# Arrange plots using arrange_plots function
arrange_plots(plot_list_Box[Set_SelectPlt], edgeOnlyXAxis = TRUE, edgeOnlyYAxis = FALSE,
              edgeOnlyXLabel = TRUE, edgeOnlyYLabel = FALSE, OneFigLegend = TRUE,
              removeTitles = TRUE, legend_position = "bottom", # legend_direction = "horizontal",
              x_tick_angle = 0,x_tick_size = 20,y_tick_size = 20,
              title_text = Name_FigTitle, removeX = TRUE,simplifyXAxis = TRUE,auto_labels = FALSE)


## Export PDF
pdf(paste0(Name_ExportFolder,"/",Name_Export,"_Sum_Box.pdf"),
    width = 22, height = 10)
arrange_plots(plot_list_Box[Set_SelectPlt], edgeOnlyXAxis = TRUE, edgeOnlyYAxis = FALSE,
              edgeOnlyXLabel = TRUE, edgeOnlyYLabel = FALSE, OneFigLegend = TRUE,
              removeTitles = TRUE, legend_position = "bottom",x_tick_angle = 0,
              x_tick_size = 28,y_tick_size = 28,
              title_text = Name_FigTitle, removeX = TRUE,simplifyXAxis = TRUE,auto_labels = FALSE)
# Bug in simplifyXAxis = TRUE
print(plot_list_Box2)
dev.off()
