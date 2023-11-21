##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

##### Set parameter* ####
Set_FindVarFea <- "vst" # "disp" # "mvp" # "vst"
Name_Title = paste0("CompDiff_mix10x5_H838_Sum500 nfeature:2000; numGeneforEst:2000; nPC:50; FindVarFea:",Set_FindVarFea)
# Name_ExportFolder = paste0("Export/#_MSINB_0625_Performance_CompDiff_H838_scKW2000_FilterAllCell_#/",Set_FindVarFea,"/","H838")
Name_ExportFolder = paste0("Export/#_MSINB_0625_Performance_CompDiff_H838_scKW2000_FilterAllCell_1119_#")

Name_Export = paste0(gsub("-", "_",Sys.Date()),"_Sum_Metics_", Set_FindVarFea)
Set_SelectPlt <- c("ASWPCA_CellType", "ARI", "Inverse_PCADepthCorr" )

##### Load data #####
Set_Path <- paste0("D:/Dropbox/##_GitHub/###_VUMC/scKWARN/scKWARN_Test/", Name_ExportFolder)
# Set_Path <- "D:/Dropbox/##_GitHub/###_VUMC/scKWARN_Test/Export/#_MSINB_0606_Performance_FindVarFea_CompDiff_#/#_Sum_mvp" # Specify directory

## Merge specific files from a folder into a Dataframe
files <- list.files(Set_Path, pattern = ".tsv", full.names = TRUE)

# # Get a list of all .tsv files in the directory that end with 'CT'.
# files <- list.files(Set_Path, pattern = "*CT*.tsv", full.names = TRUE)
files <- list.files(Set_Path, pattern = "*Sum*.tsv", full.names = TRUE)
# # ## Get a list of all files in the directory that contain "CT" and are tsv files
# # files <- list.files(Set_Path, pattern = "CT.*\\.tsv$", full.names = TRUE)

## Read each tsv file and stack them into a data frame
Sum.df <- files %>%
  map_dfr(read_tsv)

print(Sum.df) # Display the data

# Check if the column "PCADepthCorr" exists
if ("PCADepthCorr" %in% colnames(Sum.df)) {
  # If exists, create a new column "Inverse_PCADepthCorr"
  Sum.df$Inverse_PCADepthCorr <- 1 - Sum.df$PCADepthCorr
}

# ## Filter some condition
# Sum.df <- Sum.df[!Sum.df$NorMeth %in% "counts" & !Sum.df$Dataset %in% "Sum_500_H838_50", ]
Sum.df <- Sum.df[!Sum.df$Dataset %in% "Sum_500_H838_50",]
Sum.df$NorMeth <- gsub("counts", "Non-normalization", Sum.df$NorMeth)

##### Plot Box #####
library(ggplot2)
source("FUN_Plot_Box.R")
source("FUN_Plot_Arrange.R")

# Set Color
col.values <- c("Non-normalization" = "gray","scKWARN" = "blue", "RC" = "black", "scran" = "red", "SCNorm" = "#00CC00",
                "Seurat_LogNormalize" = "#159e5a","Seurat_RC" = "#2bed78","Seurat_CLR" = "#335e49",
                "sctransform" = "#B266FF", "PsiNorm" = "#FF8000")
Set_NorType <- c("Non-normalization","scKWARN", "RC", "scran",  # "SCNorm", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
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
pdf(paste0(Name_ExportFolder,"/",Name_Export,"_Sum_Box_V2.pdf"),
    width = 22, height = 10)
arrange_plots(plot_list_Box[Set_SelectPlt], edgeOnlyXAxis = TRUE, edgeOnlyYAxis = FALSE,
              edgeOnlyXLabel = TRUE, edgeOnlyYLabel = FALSE, OneFigLegend = TRUE,
              removeTitles = TRUE, legend_position = "bottom",x_tick_angle = 0,
              x_tick_size = 28,y_tick_size = 28,
              title_text = Name_FigTitle, removeX = TRUE,simplifyXAxis = TRUE,auto_labels = FALSE)
# Bug in simplifyXAxis = TRUE
print(plot_list_Box2)
dev.off()




# ######################################################################################
# ##### Plot Box (One Plot Test) #####
# source("FUN_Plot_Box.R")
#
#
# # ## Load data (Test)
# # Sum.df <- read.delim2(paste0(getwd(),"/Export/","MSINB_2023_06_05__scKW2000_Fea2000_SetCltNumTRUE_nPC50_FindVarFeadisp_V1_Sum_Test.txt"),
# #                       na.strings = "NA")
#
#
# # Replace "" with NA in Sum.df
# Sum.df <- data.frame(lapply(Sum.df, function(x) replace(x, x == "", NA)))
#
# # Set Color
# col.values <- c("counts" = "gray","scKWARN" = "blue", "RC" = "black", "scran" = "red", "SCNorm" = "#00CC00",
#                 "Seurat_LogNormalize" = "#159e5a","Seurat_RC" = "#2bed78","Seurat_CLR" = "#335e49",
#                 "sctransform" = "#B266FF", "PsiNorm" = "#FF8000")
# Set_NorType <- c("counts","scKWARN", "RC", "scran", "SCNorm",  # "SCNorm", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
#                  "sctransform", "PsiNorm")
# Set_Dataset <- c("mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53", "mix.DropSeq", "mix.10x", "mix.10x5", "mix.CELSeq")
#
# # You can then call this function like this:
# plt.box <- create_box_plot_multiXTpye(Sum.df,
#                                       x_var = "Dataset", y_var = "PCADepthCorr",
#                                       fill_var = "NorMeth", color_var = "NorMeth",
#                                       Set_fill_order = Set_NorType,
#                                       Set_x_order = Set_Dataset,
#                                       Name_XTitle = "Dataset",
#                                       AspectRatio = 0.8,
#                                       Name_YTitle = "PCA Depth Correlation",
#                                       Name_LegendTitle = "Normalization",
#                                       Set_Title = "",
#                                       color_values = col.values,
#                                       ylim_values = c(0, 1), # set y limits
#                                       show_mean = TRUE) # show mean value on boxplot
# plt.box
