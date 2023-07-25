##### Presetting ######
rm(list = ls()) # Clean variable
memory.limit(150000)

## Ref: https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html
##### Load Packages #####
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if (!requireNamespace("SeuratData", quietly = TRUE)) {install.packages("SeuratData")}; library(SeuratData)
if (!require(ComplexHeatmap)) {install.packages("ComplexHeatmap")}; library(ComplexHeatmap)
if (!require(circlize)) {install.packages("circlize")}; library(circlize)

##### Load RData #####
load(paste0(getwd(), "/Export/#_MSINB_0611_Performance_#/MSINB_2023_06_11__scKW2000_Fea2000_SetCltNumTRUE_nPC50_FindVarFeavst_V1.RData"))

seurat_list <- seurat_Nor_list

## Clean up the object
rm(list = setdiff(ls(), c("seurat_list", "summary_df", "Name_Export", "Name_ExportFolder" )))

##### Set Export & Condition #####
## Set parameters
Set_FindVarFea <- "vst" # "disp" # "mvp" # "vst" # Set variable feature option
Set_seed <- TRUE # Set seed value to TRUE

# ## Set ExportFolder
# Name_ExportFolder <- "Export"
# if (!dir.exists(Name_ExportFolder)){dir.create(Name_ExportFolder)}   ## Create new folder
#
# Name_Export <- paste0("MSINB_2023_06_11_scKW2000_Fea2000_SetCltNumTRUE_nPC50")


################################################################################
# Create an empty replicability result list
replicability <- list()
variable.features.s1 <- list()
variable.features.s2 <- list()
variable.features <- list()
seurat_split1_list <- list()
seurat_split2_list <- list()

# Your list of normalization methods
norm.methods <- c("counts","scKWARN", "RC", "scran", "sctransform", "PsiNorm", "SCNorm")

# Iterate through all the datasets and normalization methods
for (dataset in names(seurat_list)){
  for (norm in norm.methods){

    # Get the Seurat object for this dataset and normalization method
    seurat <- seurat_list[[dataset]][[norm]]
    Name_NorType <- norm

    # Replace the count slot of the Seurat object with the Normalized data
    seurat@assays[["RNA"]]@counts <- seurat@assays[["RNA"]]@data

    # Find the variable features
    seurat <- FindVariableFeatures(seurat, selection.method = Set_FindVarFea)

    # Save the variable features
    variable.features[[dataset]][[norm]] <- VariableFeatures(seurat)

    # Randomly split the Seurat object into two halves
    cell_ids <- colnames(seurat)

    if (Set_seed == TRUE) {
      set.seed(717)
    } else {
      set.seed(as.numeric(Sys.time()))  # Use the value of Sys.time() as a seed to generate different random results
    }

    split_ids <- sample(c(TRUE, FALSE), size = length(cell_ids), replace = TRUE)

    seurat.split1 <- subset(seurat, cells = cell_ids[split_ids])
    seurat.split2 <- subset(seurat, cells = cell_ids[!split_ids])


    seuratObject <- seurat.split1
    source("#_Run_Normalization_Jaccard.R")
    seurat.split1 <- seuratObject
    rm(seuratObject)

    seuratObject <- seurat.split2
    source("#_Run_Normalization_Jaccard.R")
    seurat.split2 <- seuratObject
    rm(seuratObject)

    # Skip FindVariableFeatures steps if Name_NorType equals "sctransform"
    if (Name_NorType != "sctransform") {
      ## Find the variable features for each split
      seurat.split1 <- FindVariableFeatures(seurat.split1, selection.method = Set_FindVarFea) # nfeatures = 2000
      seurat.split2 <- FindVariableFeatures(seurat.split2, selection.method = Set_FindVarFea) # nfeatures = 2000

    }


    # Extract the variable features
    variable.features.split1 <- VariableFeatures(seurat.split1)[1:500] %>% na.omit()
    variable.features.split2 <- VariableFeatures(seurat.split2)[1:500] %>% na.omit()

    variable.features.s1[[dataset]][[norm]] <- variable.features.split1
    variable.features.s2[[dataset]][[norm]] <- variable.features.split2

    # Calculate the replicability
    intersection <- length(intersect(variable.features.split1, variable.features.split2))
    union <- length(union(variable.features.split1, variable.features.split2))
    replicability[[dataset]][[norm]] <- intersection / union

    # Record the result
    seurat_split1_list[[dataset]][[norm]] <- seurat.split1
    seurat_split2_list[[dataset]][[norm]] <- seurat.split2

    ## Clean up the obj
    rm(seurat, seurat.split1, seurat.split2, variable.features.split1, variable.features.split2)
  }
}


# Initialize the reproducibility result list
reproducibility <- list()

# Iterate through all the normalization methods
for (norm in norm.methods){

  # Initialize the list for this normalization method
  reproducibility[[norm]] <- list()

  # Iterate through all the combinations of datasets
  for (i in 1:(length(names(seurat_list))-1)){
    for (j in (i+1):length(names(seurat_list))){

      # Get the variable features for these two datasets and this normalization method
      vf1 <- variable.features[[names(seurat_list)[i]]][[norm]]
      vf2 <- variable.features[[names(seurat_list)[j]]][[norm]]

      # Calculate the reproducibility
      intersection <- length(intersect(vf1, vf2))
      union <- length(union(vf1, vf2))
      reproducibility[[norm]][[paste(names(seurat_list)[i], names(seurat_list)[j], sep="-")]] <- intersection / union
    }
  }
}

# Combine the results into a data frame
replicability_df <- do.call(rbind, lapply(replicability, function(x) setNames(data.frame(t(unlist(x))), names(x))))
# view(replicability_df)
reproducibility_df <- do.call(rbind, lapply(reproducibility, function(x) setNames(data.frame(t(unlist(x))), names(x))))
reproducibility_df <- reproducibility_df %>% t()
# view(reproducibility_df)


################################################################################
##### Plot Heatmap #####
replicability_df <- as.matrix(replicability_df)
reproducibility_df <- as.matrix(reproducibility_df)
# Create color mapping function
col_fun <- colorRamp2(c(0, 1), c("lightblue", "darkblue"))

replicability_heatmap <- Heatmap(replicability_df,
                                 name = "replicability",
                                 col = col_fun,
                                 cell_fun = function(j, i, x, y, width, height, fill) {
                                   grid.text(sprintf("%.3f", replicability_df[i, j]), x, y, gp = gpar(fontsize = 10))
                                 })

reproducibility_heatmap <- Heatmap(reproducibility_df,
                                   name = "reproducibility",
                                   col = col_fun,
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     grid.text(sprintf("%.3f", reproducibility_df[i, j]), x, y, gp = gpar(fontsize = 10))
                                   })

# Draw heatmaps stacked vertically
draw(replicability_heatmap %v% reproducibility_heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")

draw(replicability_heatmap, heatmap_legend_side = "bot")
draw(reproducibility_heatmap, heatmap_legend_side = "bot")

#### Combine Data ####
# Combine the matrices vertically
all_data <- rbind(replicability_df, reproducibility_df)

# Create color mapping function
col_fun <- colorRamp2(c(0, 1), c("lightblue", "darkblue"))

# Create the combined heatmap
all_heatmap <- Heatmap(all_data,
                       name = "Jaccard index",
                       col = col_fun,
                       show_row_names = TRUE,
                       show_row_dend = FALSE,
                       cluster_rows = TRUE, # FALSE,
                       row_split = factor(rep(c("replicability", "reproducibility"),
                                              times = c(nrow(replicability_df), nrow(reproducibility_df)))),
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.3f", all_data[i, j]), x, y, gp = gpar(fontsize = 10))
                       })

# Draw the heatmap with legend on the right and a title
#(Error)# draw(all_heatmap, heatmap_legend_side = "right", main = "Your title")
draw(all_heatmap,column_title = paste0("FindVariableFeatures: ",Set_FindVarFea),
     column_title_gp = gpar(fontsize = 16))


## Export PDF
pdf(paste0(Name_ExportFolder,"/","#_",Name_Export,"_JaccardIndex_FindVarFea",Set_FindVarFea,"_FixSeed",Set_seed,".pdf"),
    width = 12, height = 7)
draw(all_heatmap,column_title = paste0("FindVariableFeatures: ",Set_FindVarFea),
     column_title_gp = gpar(fontsize = 16))
dev.off()

#################################################################################
## Summarize
library(magrittr)

# Define the group
group1 <- c("mix.CELSeq51", "mix.CELSeq52", "mix.CELSeq53")

# Mean_all & Mean_grouped
replicability_Sum_df <- replicability_df %>%
  {tmp_df <- .; rbind(tmp_df, mean_all = colMeans(tmp_df))} %>%
  {tmp_df <- .; rbind(tmp_df, mean_group1 = colMeans(tmp_df[group1, ]))} %>%
  {tmp_df <- .; rbind(tmp_df, mean_grouped = colMeans(tmp_df[setdiff(rownames(tmp_df), c(group1, "mean_all")), ]))}

#################################################################################
#### Export #####
## Export result df
write.table(data.frame(Jaccard = row.names(replicability_Sum_df), replicability_Sum_df),
            file=paste0(Name_ExportFolder,"/",Name_Export,"_Replicability_Sum.tsv"),
            quote = FALSE,row.names = FALSE,col.names = TRUE, na = "",sep = '\t')


## Export RData
save.image(paste0(Name_ExportFolder,"/",Name_Export,"_JaccardIndex_FindVarFea",Set_FindVarFea,"_FixSeed",Set_seed,".RData"))


#################################################################################
