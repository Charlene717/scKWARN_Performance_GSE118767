# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)

## Record time set
Rec_Time_Point.lt <- list()
Rec_Time_Spend.lt <- list()

Rec_Time_Point.lt[["Start_Time"]] <- Sys.time() # %>% as.character()


##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if (!requireNamespace("SeuratData", quietly = TRUE)) {install.packages("SeuratData")}; library(SeuratData)

# ## Call function
# source("#_FUN_TimeRec.R")

Rec_Time_Point.lt[["Load_Packages"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Load_Packages"]] <- Rec_Time_Point.lt[["Load_Packages"]]-Rec_Time_Point.lt[["Start_Time"]]

##### Set Export & Condition #####
## Set parameters
# Set_numGeneforEst = 2000
Set_nfeatures = 2000
Num_PCA = 50
Set_TarClusterNum <- TRUE

Set_logPlusOne = TRUE
Set_logDepth = TRUE
# Set_Filter_Cell = TRUE

# Set_FindVarFea <- "vst" # "disp" # "mvp" # "vst"

## Set Normalization method
# Set_NorType <- c("counts","scKWARN", "RC", "scran", "SCNorm", # "SCNorm", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
#                  "sctransform", "PsiNorm")
Set_sctransform_Res = TRUE

## Set the metrics to be included
Set_Metrics <- c("PCADepthCorr", "ASWPCA_CellType","ARI") # Set_Metrics <- c("PCADepthCorr", "ASWPCA_CellType","ASWPCA_Cluster","ARI", "Purity", "ROGUE_CellType", "ROGUE_Cluster")

## Set export integrated plot
Set_SelectPlt_Multi <- c("ASWPCA_CellType", "Inverse_PCADepthCorr", "ARI" )
Set_Dataset_group1 <- c("mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53") # Define Dataset group


## Set names
Name_CP <- "MSINB"
# Set_Test_Type <- "CompDiff" # "CompDiff" or ""
# Set_Seed <- 123
# Name_Test <- "V1"

if(Set_Test_Type == "CompDiff"){
  # Set_TotalCellNum <- 800
  Set_YLimCellNum <- Set_TotalCellNum
  # Set_Main_CellLine = "A549"
  # Set_Main_FltCellLine = "ALL" # Set_Main_FltCellLine = Set_Main_CellLine
  Name_DataSet <- "mix10x5"

  Name_Sup <- paste0("CompDiff_",Name_DataSet,"_",Set_Main_CellLine,"_Sum",Set_TotalCellNum) # Name_Sup <- "Test" # Name_Sup <- "CompDiff"
  Set_Ratio <- c(0.9, 0.8, 0.6, 0.4, 0.2, 0.1)

}else{
  Set_YLimCellNum <- 1400
  # Define the datasets to include in the new seurat_list
  Set_Dataset <- c( "mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53", "mix.DropSeq", "mix.10x",
                    "mix.10x5", "mix.CELSeq" )

  Name_Sup <- ""
}

Name_FigTitle <- paste0(Name_Sup,"  ","nfeatures:",Set_nfeatures,"; numGeneforEst:", Set_numGeneforEst, # "; Log(Nor+1):",Set_logPlusOne,
                        "; nPC:",Num_PCA,"; FindVarFea:",Set_FindVarFea,"; Res:", Set_sctransform_Res)


Name_Note <- paste0(Name_Sup, "_scKW",Set_numGeneforEst, #"_LogNp1",Set_logPlusOne,
                    "_Fea",Set_nfeatures, "_SetCltNum",Set_TarClusterNum,
                    "_nPC",Num_PCA,"_FindVarFea",Set_FindVarFea,"_",Name_Test)

Name_Export <- paste0(Name_CP,"_",
                     gsub("-", "_",Sys.Date()),
                     "_",Name_Note) # as.numeric(Sys.time(), units = "secs") %>% trunc())

## Set ExportFolder
Name_ExportFolder <- "Export"
if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder

##### Load data* #####
seurat_list <- list()
source("Dataset_PsiNorm.R")

## Add "Cell_Type" col to metadata in Seurat object
source("Run_Rename_Seurat_MetadataColname.R")


if(Set_Test_Type == "CompDiff"){
  if (Set_Filter_Cell) {
    ## Filter Cell
    source("Dataset_PsiNorm_mix10x5_Filter.R")
    seurat_list[["mix.10x5"]] <- seurat_object_filtered
    ## Case of Composition Differences
    source("Run_Rename_Seurat_MetadataColname.R")
  }
  source("##_Create_Compositional_Differences.R")

}else{
  # Remove datasets not in Set_Dataset from the seurat_list
  seurat_list <- seurat_list[names(seurat_list) %in% Set_Dataset]

}

#### QC ####
source("Run_Seurat_QC.R")

Rec_Time_Point.lt[["Load_data"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Load_data"]] <- Rec_Time_Point.lt[["Load_data"]] - Rec_Time_Point.lt[["Load_Packages"]]

#### Normalizing the data ####
source("#_FUN_NorMeth.R") # Call function
source("#_Run_Normalization_SeuratList.R")

Rec_Time_Point.lt[["Normalization"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Normalization"]] <- Rec_Time_Point.lt[["Normalization"]] - Rec_Time_Point.lt[["Load_data"]]

##### Clustering & DR #####
source("FUN_Seurat_Set_ClusterNum.R")
# Initialize the result list
seurat_Result_list <- list()

# Loop through each main element in seurat_Nor_list (e.g., mix.CELSeq51, mix.CELSeq52)
for (main_name in names(seurat_Nor_list)) {
  # Initialize a list to hold the processed Seurat objects for this main element
  main_result_list <- list()

  # Loop through each sub element (e.g., RC, scran, sctransform, PsiNorm, scKWARN)
  for (sub_name in names(seurat_Nor_list[[main_name]])) {
    try({
    # Get the Seurat object
    seuratObject <- seurat_Nor_list[[main_name]][[sub_name]]
    # Replace the count slot of the Seurat object with the Normalized data
    seuratObject@assays[["RNA"]]@counts <- seuratObject@assays[["RNA"]]@data

    # Apply the preprocessing steps
    # seuratObject <- NormalizeData(seuratObject) # Different normalization methods have been tried earlier

    # Skip FindVariableFeatures and ScaleData steps if sub_name equals "sctransform"
    if (sub_name != "sctransform") {
      seuratObject <- FindVariableFeatures(seuratObject,nfeatures = Set_nfeatures,
                                           selection.method = Set_FindVarFea) %>%
                      ScaleData()
    }

    seuratObject <- RunPCA(seuratObject) %>% # Dimensionality reduction
                    FindNeighbors(dims = 1:Num_PCA)

    # Calculate target_clusters if needed
    if (exists("Set_TarClusterNum") && Set_TarClusterNum) { # Use the function to perform clustering with a target number of clusters
      seuratObject <- cluster_Seurat_to_target(seuratObject, target_clusters = length(unique(seuratObject@meta.data$Cell_Type))) # count unique cell types in metadata
    } else { # Perform standard Seurat clustering (Cluster cells based on PCA)
      seuratObject <- FindClusters(seuratObject)
    }

    # Dimensionality reduction for visualization
    seuratObject <- RunUMAP(seuratObject, dims = 1:Num_PCA) %>%
                    RunTSNE(dims = 1:Num_PCA)

    # Save the processed Seurat object back into the main_result_list
    main_result_list[[sub_name]] <- seuratObject
    })
  }

  # Save the main_result_list into the seurat_Result_list
  seurat_Result_list[[main_name]] <- main_result_list
}

# rm(seurat_Nor_list)
Rec_Time_Point.lt[["Clt_DR"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Clt_DR"]] <- Rec_Time_Point.lt[["Clt_DR"]] - Rec_Time_Point.lt[["Normalization"]]

##### Metrics #####
## Calculate Metrics
library(cluster)
source("#_FUN_Metric.R")

# # Set the metrics to be included
# Set_Metrics <- c("PCADepthCorr","ARI", "ASWPCA_CellType", "ASWPCA_Cluster",
#                  "Purity", "ROGUE_CellType", "ROGUE_Cluster")

# Iterate over all datasets in the list
for (dataset_name in names(seurat_Result_list)) {
  for (method_name in names(seurat_Result_list[[dataset_name]])) {
    # Get the current Seurat object
    seuratObject <- seurat_Result_list[[dataset_name]][[method_name]]

    # Iterate over the metrics
    for (metric in Set_Metrics) {

      if (metric == "PCADepthCorr") {
        numPCs <- 2  # Set this variable to determine the number of PCs
        PCADepthCorr <- calculate_PCADepthCorr(seuratObject, numPCs, logDepth = Set_logDepth)
        seuratObject@misc[["PCADepthCorr"]] <- PCADepthCorr

      } else if (metric == "ASWPCA_CellType") {
        results_ASWPCA_CellType <- calculate_silhouette(seurat_object = seuratObject, reduction_method = "pca", distance_method = "euclidean", cluster_method = "Cell_Type", pca_components = Num_PCA)
        seuratObject@misc[["ASWPCA_CellType"]] <- results_ASWPCA_CellType
      } else if (metric == "ASWPCA_Cluster") {
        results_ASWPCA_Cluster <- calculate_silhouette(seurat_object = seuratObject, reduction_method = "pca", distance_method = "euclidean", cluster_method = "seurat_clusters", pca_components = Num_PCA)
        seuratObject@misc[["ASWPCA_Cluster"]] <- results_ASWPCA_Cluster

      } else if (metric == "Purity") {
        cluster_purity.df <- calculate_purity(seurat_object = seuratObject, annotation_name = "Cell_Type", cluster_name = "seurat_clusters")
        seuratObject@misc[["Purity"]] <- cluster_purity.df

      } else if (metric == "ROGUE_CellType") {
        try({
          cluster_ROGUE.lt <- calculate_rogue_value(seuratObject, labels = "Cell_Type", samples = "Cell_Type")
          seuratObject@misc[["ROGUE_CellType"]] <- cluster_ROGUE.lt
        })
      } else if (metric == "ROGUE_Cluster") {
        try({
          celltype_ROGUE.lt <- calculate_rogue_value(seuratObject, labels = "seurat_clusters", samples = "seurat_clusters")
          seuratObject@misc[["ROGUE_Cluster"]] <- celltype_ROGUE.lt
        })

      } else if (metric == "ARI") {
        try({
          seuratObject@misc[["ARI"]] <- calculate_ari(seuratObject, "Cell_Type", "seurat_clusters")
        })
      }

    }

    # Save the Seurat object back into the list
    seurat_Result_list[[dataset_name]][[method_name]] <- seuratObject
  }
}

## Plot CellCount list
source("FUN_create_plots_CellCount_list.R")

# Use the function to create a list of plots with reverse order and sum of cells in the title
plots_CellCount_list <- create_plots_CellCount_list(seurat_list, reverse_order = TRUE, add_sum_cells = TRUE, YLimCellNum = Set_YLimCellNum)

# Combine all plots into a grid
combined_plot_CellCount <- grid.arrange(grobs = plots_CellCount_list, ncol = 3)  # Adjust the number of columns if necessary

## Export PDF
pdf(paste0(Name_ExportFolder,"/",Name_Export,"_CellCount.pdf"),
    width = 12, height = 10)
combined_plot_CellCount <- grid.arrange(grobs = plots_CellCount_list, ncol = 3)  # Adjust the number of columns if necessary
dev.off()


## Export TIFF
tiff(paste0(Name_ExportFolder,"/",Name_Export,"_CellCount.tiff"),
     width = 800, height = 700)
combined_plot_CellCount <- grid.arrange(grobs = plots_CellCount_list, ncol = 3)  # Adjust the number of columns if necessary
dev.off()

rm(PCADepthCorr, cluster_purity.df, results_ASWPCA_Cluster,results_ASWPCA_CellType,
   cluster_ROGUE.lt,celltype_ROGUE.lt)

Rec_Time_Point.lt[["Metrics"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Metrics"]] <- Rec_Time_Point.lt[["Metrics"]] - Rec_Time_Point.lt[["Clt_DR"]]

##### Summarize the result #####
#### Summarize all Value ####
# Initialize an empty data frame
summary_df <- data.frame(Dataset = character(),NorMeth = character(),stringsAsFactors = FALSE)

# Iterate over all datasets in the list
for (dataset in names(seurat_Result_list)) {
  # Iterate over all normalization methods in the current dataset
  for (normethod in names(seurat_Result_list[[dataset]])) {
    # Get the current Seurat object
    seuratObject <- seurat_Result_list[[dataset]][[normethod]]

    # Create a new row with dataset and normalization method
    new_row <- data.frame(Dataset = dataset, NorMeth = normethod, stringsAsFactors = FALSE)

    # Check and add the desired metrics
    for (metric in Set_Metrics) {
      metric_value <- NA  # Initialize the metric value to NA

      if (metric %in% names(seuratObject@misc)) {
        if (metric == "PCADepthCorr") {
          metric_value <- seuratObject@misc[[metric]][["max_correlation"]]
        } else if (metric %in% c("ASWPCA_CellType", "ASWPCA_Cluster")) {
          metric_value <- seuratObject@misc[[metric]][["average_silhouette_width"]]
        } else if (metric == "Purity") {
          metric_value <- mean(seuratObject@misc[[metric]][["Purity"]])
        } else if (metric %in% c("ROGUE_CellType","ROGUE_Cluster")) {
          metric_value <- seuratObject@misc[[metric]][["avg_rogue.value"]]
        } else if (metric == "ARI") {
          metric_value <- seuratObject@misc[[metric]]
        }
      }

      new_row[[metric]] <- metric_value
    }

    # Get the number of genes and cells
    new_row$GeneNum <- length(rownames(seuratObject@assays$RNA@counts))
    new_row$CellNum <- length(colnames(seuratObject@assays$RNA@counts))

    if (exists("Set_TarClusterNum") && Set_TarClusterNum) {
      # Get the number and resolution of cluster
      new_row$clustersNum <- seuratObject@misc[["Rec_clustersNum"]]$num_clusters
      new_row$clustersResolution <- seuratObject@misc[["Rec_clustersNum"]]$resolution
      new_row$TarclustersNum <- seuratObject@misc[["Rec_clustersNum"]]$target_clustersNum
      new_row$Record <- Name_Export

    }else{
      new_row$clustersNum <- length(unique(seuratObject@meta.data$seurat_clusters))
      new_row$TarclustersNum <- length(unique(seuratObject@meta.data$Cell_Type))
    }


    # Append the new row to the summary data frame
    summary_df <- rbind(summary_df, new_row)
  }
}

# Reset NA values to "NA" for character columns
summary_df[] <- lapply(summary_df, function(x) if (is.character(x)) replace(x, is.na(x), "NA") else x)
# summary_df$Dataset <- gsub("_", ".", summary_df$Dataset)

# Check if the column "PCADepthCorr" exists
if ("PCADepthCorr" %in% colnames(summary_df)) {
  # If exists, create a new column "Inverse_PCADepthCorr"
  summary_df$Inverse_PCADepthCorr <- 1 - summary_df$PCADepthCorr
  Set_Metrics <- c(Set_Metrics,"Inverse_PCADepthCorr")
}


#### Summarize of Mean Value ####
library(tidyverse)
# Set_Dataset_group1 <- c("mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53") # Define Dataset group

# Convert the data into a format that is easier to work with
tidy_df <- summary_df %>%
  pivot_longer(cols = Set_Metrics, names_to = "Metric", values_to = "Value") # cols = c(PCADepthCorr, ASWPCA_CellType, ASWPCA_Cluster, ARI, Purity, ROGUE_CellType, ROGUE_Cluster),

## Compute the average of all Datasets
summary_mean_all_df <- tidy_df %>%
  group_by(NorMeth, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = NorMeth, values_from = Value)

## Compute average of grouped Datasets
# First, calculate the mean of Set_Dataset_group1
group1_tidy_df <- tidy_df %>%
  filter(Dataset %in% Set_Dataset_group1) %>%
  group_by(NorMeth, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  mutate(Dataset = "group1")

# Then, get other Dataset
other_tidy_df <- tidy_df %>%
  filter(!Dataset %in% Set_Dataset_group1)

# Merge two DataFrames
grouped_tidy_df <- bind_rows(group1_tidy_df, other_tidy_df)

# Compute mean of grouped Datasets
summary_mean_grouped_df <- grouped_tidy_df %>%
  group_by(NorMeth, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  pivot_wider(names_from = NorMeth, values_from = Value)

Rec_Time_Point.lt[["SumResult"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["SumResult"]] <- Rec_Time_Point.lt[["SumResult"]] - Rec_Time_Point.lt[["Metrics"]]

##### Visualization #####
## Basic setting
library(ggplot2)
# Set Color
col.values <- c("counts" = "gray","scKWARN" = "blue", "RC" = "black", "scran" = "red", "SCNorm" = "#00CC00",
                "Seurat_LogNormalize" = "#159e5a","Seurat_RC" = "#2bed78","Seurat_CLR" = "#335e49",
                "sctransform" = "#B266FF","PsiNorm" = "#FF8000")

# Define fixed parameters for plot_bubble function
fixed_params <- list(x = "Dataset", fill = "NorMeth", col.values = col.values,
                     x_order = Set_Dataset, fill_order = Set_NorType,
                     NameX = "Dataset", NameFill = "Normalization")

# Define list for variable parameters with names
variable_params_list <- list(
  ASWPCA_CellType = list(df = summary_df[!grepl("pbmc3k", summary_df$Dataset), ], y = "ASWPCA_CellType",
                         NameY = "ASW (Cell Type)", Set_TitleName = Name_FigTitle),
  ASWPCA_Cluster = list(df = summary_df, y = "ASWPCA_Cluster",
                        NameY = "ASW (Cluster)", Set_TitleName = Name_FigTitle),
  Purity = list(df = summary_df[!grepl("pbmc3k", summary_df$Dataset), ], y = "Purity",
                NameY = "Purity", Set_TitleName = Name_FigTitle),
  PCADepthCorr = list(df = summary_df, y = "PCADepthCorr",
                      NameY = "Correlation", Set_TitleName = paste0(Name_FigTitle,"; Log(Depth):",Set_logDepth)),
  Inverse_PCADepthCorr = list(df = summary_df, y = "Inverse_PCADepthCorr",
                              NameY = "1-Correlation", Set_TitleName = paste0(Name_FigTitle,"; Log(Depth):",Set_logDepth)),
  ROGUE_CellType = list(df = summary_df[!grepl("pbmc3k", summary_df$Dataset), ], y = "ROGUE_CellType",
                        NameY = "ROGUE (Cell Type)", Set_TitleName = Name_FigTitle),
  ROGUE_Cluster = list(df = summary_df, y = "ROGUE_Cluster",
                       NameY = "ROGUE (Cluster)", Set_TitleName = Name_FigTitle),
  ARI = list(df = summary_df[!grepl("pbmc3k", summary_df$Dataset), ], y = "ARI",
             NameY = "ARI", Set_TitleName = Name_FigTitle)
)

# Set_SelectPlt_Multi <- c("ASWPCA_CellType","Purity","PCADepthCorr","ROGUE_CellType" , "ARI" )

variable_params_list <- variable_params_list[Set_Metrics]


#### Bar plot ####
# Use lapply to apply plot_bar function
source("FUN_Plot_Bar.R")
source("FUN_Plot_Arrange.R")
plot_list_Bar <- lapply(variable_params_list, function(params) {
  args <- modifyList(fixed_params, params)
  do.call(plot_bar, args)
})
print(plot_list_Bar)

## Export PDF
pdf(paste0(Name_ExportFolder,"/",Name_Export,"_Sum_Bar.pdf"),
    width = 12, height = 10)
arrange_plots(plot_list_Bar[Set_SelectPlt_Multi], edgeOnlyXAxis = TRUE, edgeOnlyYAxis = FALSE,
              edgeOnlyXLabel = TRUE, edgeOnlyYLabel = FALSE, OneFigLegend = TRUE,
              removeTitles = TRUE, legend_position = "bottom",x_tick_angle = 45,x_tick_size = 12,
              title_text = Name_FigTitle,removeX = TRUE,simplifyXAxis = TRUE,auto_labels = FALSE) %>% print()

print(plot_list_Bar)
dev.off()

## Export TIFF
# Loop through each plot and export it
lapply(names(plot_list_Bar), function(name) {
  tiff(paste0(Name_ExportFolder, "/", Name_Export, "_Bar_", name, ".tiff"),
       width = 800, height = 700)
  print(plot_list_Bar[[name]])
  dev.off()
})


##### Export #####
## Export result df
write.table(summary_df,
            file=paste0(Name_ExportFolder,"/",Name_Export,"_Sum.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

## Export summary_mean_all_df
write.table(summary_mean_all_df,
            file=paste0(Name_ExportFolder,"/",Name_Export,"_Sum_MeanAll.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')

## Export summary_mean_grouped_df
write.table(summary_mean_grouped_df,
            file=paste0(Name_ExportFolder,"/",Name_Export,"_Sum_MeanGrouped.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')


Rec_Time_Point.lt[["End_Time"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["Sum"]] <- Rec_Time_Point.lt[["End_Time"]] - Rec_Time_Point.lt[["Start_Time"]]


## Record version and sessionInfo
info_output <- c("##_R Version Information:", capture.output(version), "",
                 "##_Session Information:", capture.output(sessionInfo()))

writeLines(info_output, paste0(Name_ExportFolder,"/",Name_Export,"_Version_and_Session_Info_Sum.txt"))


# ## Export RData
# save.image(paste0(Name_ExportFolder,"/",Name_Export,".RData"))

Rec_Time_Point.lt[["SaveRData"]] <- Sys.time() # %>% as.character()
Rec_Time_Spend.lt[["SaveRData"]] <- Rec_Time_Point.lt[["SaveRData"]] - Rec_Time_Point.lt[["End_Time"]]


