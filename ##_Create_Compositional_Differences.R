## Ref: https://satijalab.org/seurat/articles/integration_introduction.html

# ##### Presetting ######
# rm(list = ls()) # Clean variable
# memory.limit(150000)
# Set_Seed <- 717

##### Load required libraries #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)


##### Load data #####
source("Run_Dataset_Selection_and_Integration.R")


##### Recompose_seurat #####
recompose_seurat <- function(seurat_object, cell_type_field, sampling_mode = "equal", # sampling_mode = "equal" or "random"
                             set_comp) {

  # Cell type counts
  cell_types <- seurat_object@meta.data[[cell_type_field]]
  cell_type_counts <- table(cell_types)

  # Check if the mode is "fixed total cell count"
  fixed_total_count_mode <- names(set_comp)[1] == "Sum"

  if (fixed_total_count_mode) {
    # Compute the number of cells to sample for each cell type
    cell_type_samples <- set_comp[-1] * set_comp[1]
    # Compute remaining number of cells to sample
    remaining_cells <- set_comp[1] - sum(cell_type_samples)
  } else {
    # Use the numbers in set_comp directly
    cell_type_samples <- set_comp
    # No remaining cells to sample
    remaining_cells <- 0
  }

  # Verify if all specified cell types exist in the data
  specified_cell_types <- names(cell_type_samples)
  missing_cell_types <- setdiff(specified_cell_types, names(cell_type_counts))
  if (length(missing_cell_types) > 0) {
    stop(paste("The following cell types were not found in the data:", paste(missing_cell_types, collapse = ", ")))
  }

  # Check if sampling exceeds cell type counts
  for (cell_type in names(cell_type_samples)) {
    if (cell_type_samples[cell_type] > cell_type_counts[cell_type]) {
      stop("Sampling count exceeds the number of cells available for cell type:", cell_type)
    }
  }

  if (remaining_cells > 0 && remaining_cells > sum(cell_type_counts[names(cell_type_samples)])) {
    stop("Sampling count exceeds the number of cells available for remaining cell types.")
  }

  # Create a container for the sampled cells
  sampled_cells <- list()

  for (cell_type in names(cell_type_samples)) {
    # Sample cells
    cells_to_sample <- which(cell_types == cell_type)
    if(length(cells_to_sample) < cell_type_samples[cell_type]) {
      stop("Not enough samples for cell type:", cell_type)
    }
    set.seed(ifelse(exists("Set_Seed"), Set_Seed, as.integer(Sys.time())))
    sampled_cells[[cell_type]] <- sample(cells_to_sample, cell_type_samples[cell_type], replace = FALSE)
  }

  # Sample remaining cells from the cell types not specified in set_comp
  remaining_cell_types <- setdiff(names(cell_type_counts), names(cell_type_samples))

  if (remaining_cells > 0) {
    remaining_cells_to_sample <- which(cell_types %in% remaining_cell_types)

    if (sampling_mode == "random") {
      # Random sampling
      if(length(remaining_cells_to_sample) < remaining_cells) {
        stop("Not enough samples for remaining cell types.")
      }
      set.seed(ifelse(exists("Set_Seed"), Set_Seed, as.integer(Sys.time())))
      sampled_cells[["remaining"]] <- sample(remaining_cells_to_sample, remaining_cells, replace = FALSE)

    } else if (sampling_mode == "equal") {
      # Equal distribution among remaining cell types
      cells_per_type <- floor(remaining_cells / length(remaining_cell_types))
      extra_cells <- remaining_cells - (cells_per_type * length(remaining_cell_types))

      for (cell_type in remaining_cell_types) {
        cells_to_sample <- which(cell_types == cell_type)
        num_cells_to_sample <- if (extra_cells > 0) cells_per_type + 1 else cells_per_type
        extra_cells <- extra_cells - 1

        if(length(cells_to_sample) < num_cells_to_sample) {
          stop("Not enough samples for cell type:", cell_type)
        }
        set.seed(ifelse(exists("Set_Seed"), Set_Seed, as.integer(Sys.time())))
        sampled_cells[[cell_type]] <- sample(cells_to_sample, num_cells_to_sample, replace = FALSE)
      }

    } else {
      stop(paste("Unknown sampling_mode:", sampling_mode))
    }
  }

  # Combine all sampled cells
  all_sampled_cells <- unlist(sampled_cells)

  # Create new Seurat object with sampled cells
  new_seurat_object <- subset(seurat_object, cells = all_sampled_cells)

  return(new_seurat_object)
}

# ## Test function
# seuratObject <- seurat_list[["mix.10x5"]]
# TestRecom1 <- recompose_seurat(seuratObject, cell_type_field = "Cell_Type",
#                                set_comp = c(Sum=500,A549 = 0.5,H1975 = 0.1))
# calculate_counts(TestRecom1, field_names = "Cell_Type")
#
#
# TestRecom2 <- recompose_seurat(seuratObject, cell_type_field = "Cell_Type",
#                                set_comp = c(A549 = 200,H1975 = 100))
# calculate_counts(TestRecom2, field_names = "Cell_Type")
#
# TestRecom3 <- recompose_seurat(seuratObject, cell_type_field = "Cell_Type",
#                                set_comp = c(A549 = 500,H1975 = 100))
# calculate_counts(TestRecom3, field_names = "Cell_Type")


#### Create New seurat_list #####
# Set_TotalCellNum <- 1000
# Set_Main_CellLine = "A549"
# Set_Ratio <- c(0.9, 0.8, 0.6, 0.4, 0.2, 0.1)

seurat_list <- lapply(Set_Ratio, function(ratio) {
  # Create the composition vector for each ratio
  set_comp_values <- c("Sum" = Set_TotalCellNum)
  set_comp_values[Set_Main_CellLine] = ratio

  # Create the Seurat object using the modified recompose_seurat function
  new_seurat_object <- recompose_seurat(seuratObject,
                                        cell_type_field = "Cell_Type",
                                        set_comp = set_comp_values)

  return(new_seurat_object)
})


# Assign names to the Seurat objects in the list
names(seurat_list) <- paste("Sum", Set_TotalCellNum, Set_Main_CellLine, round(Set_Ratio*Set_TotalCellNum), sep = "_")
Set_Dataset <- paste("Sum", Set_TotalCellNum, Set_Main_CellLine, round(rev(Set_Ratio)*Set_TotalCellNum), sep = "_")



# #### Check Results ####
# ## Plot CellCount list
# source("FUN_create_plots_CellCount_list.R")
#
# # Use the function to create a list of plots with reverse order and sum of cells in the title
# plots_CellCount_list <- create_plots_CellCount_list(seurat_list, reverse_order = TRUE, add_sum_cells = FALSE, YLimCellNum = Set_YLimCellNum)
#
# # Combine all plots into a grid
# combined_plot_CellCount <- grid.arrange(grobs = plots_CellCount_list, ncol = 3)  # Adjust the number of columns if necessary
#
# ## Export
# # Check if Name_ExportFolder and Name_Export exist, if not assign default values
# if (!exists("Name_ExportFolder")) {
#   Name_ExportFolder <- "Export_Temp"
#   if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
#
# }
#
# if (!exists("Name_Export")) {
#   Name_Export <- gsub("-", "_", Sys.Date())
# }
#
# ## Export PDF
# pdf(paste0(Name_ExportFolder,"/",Name_Export,"_CompDiffCheck.pdf"),
#     width = 11, height = 8)
#
# combined_plot_CellCount <- grid.arrange(grobs = plots_CellCount_list, ncol = 3)  # Adjust the number of columns if necessary
#
# dev.off()
#
# tiff(paste0(Name_ExportFolder, "/", Name_Export, "_CompDiffCheck.tiff"),
#      width = 800, height = 700)
# combined_plot_CellCount <- grid.arrange(grobs = plots_CellCount_list, ncol = 3)  # Adjust the number of columns if necessary
#
# dev.off()
#
# rm(plots_CellCount_list, combined_plot_CellCount)

################################################################################
################################# Back up code #################################
# #### Count cell number ####
# source("#_FUN_Metric.R")
# # Num_Cell.df <- calculate_counts(seuratObject, field_names = "Cell_Type")
# # Num_Cell.df
#
# seurat_list <- lapply(seurat_list, function(seurat_object) {
#   # Calculate counts
#   counts <- calculate_counts(seurat_object, field_names = "Cell_Type")
#
#   # Store the result in the misc slot
#   seurat_object@misc[["CountCell"]] <- counts
#
#   return(seurat_object)
# })


#### Create_Compositional_Differences ####
# ## Back up ##
#
# ## Set Multiple cell type
# Set_TotalCellNum <- 100
#
# # A list containing named vectors. The names of each vector are the cell lines and values are the ratios
# Set_CellLine_Ratios <- list(
#   Scenario1 = c("H1975" = 0.5, "H2228" = 0.3, "HCC827" = 0.2),
#   Scenario2 = c("H1975" = 0.4, "H2228" = 0.4, "HCC827" = 0.2),
#   Scenario3 = c("H1975" = 0.3, "H2228" = 0.4, "HCC827" = 0.3)
# )
#
# source("#_FUN_Metric.R")
# Num_Cell.df <- calculate_counts(seuratObject, field_names = "Cell_Type")
# Num_Cell.df
#
# seurat_list <- list()
# # Loop through the Set_CellLine_Ratios list
# for (scenario in names(Set_CellLine_Ratios)) {
#   cell_line_ratios <- Set_CellLine_Ratios[[scenario]]
#
#   # Create the composition vector for each scenario
#   set_comp_values <- c("Sum" = Set_TotalCellNum)
#   set_comp_values <- c(set_comp_values, cell_line_ratios)
#
#   # Create the Seurat object using the modified recompose_seurat function
#   new_seurat_object <- recompose_seurat(seuratObject,
#                                         cell_type_field = "Cell_Type",
#                                         set_comp = set_comp_values)
#
#   # Add the new Seurat object to the list
#   seurat_list[[scenario]] <- new_seurat_object
# }
#
# # The names of the Seurat objects in the list are already assigned by the names of Set_CellLine_Ratios
