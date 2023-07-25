## Unify the name of the Cell Type in the column of Metadata in the Seurat object

## Add "Cell_Type" col to metadata in Seurat object
# Array of possible column names
possible_cols <- c("cell_line_demuxlet", "seurat_annotations")

# Regular expression to match columns with "cell" and "type" (case-insensitive)
regex <- "(?i)cell.*type"

# Iterate over each Seurat object in the list
for(i in seq_along(seurat_list)) {
  # Get the current Seurat object
  seuratObject <- seurat_list[[i]]

  # Check each possible column
  for(col in possible_cols) {
    if(col %in% colnames(seuratObject@meta.data)) {
      # If it exists, copy the values to the new "Cell_Type" field
      seuratObject@meta.data$Cell_Type <- seuratObject@meta.data[, col]
    }
  }

  # Check columns matching the regex
  matching_cols <- grep(regex, colnames(seuratObject@meta.data), value = TRUE)
  for(col in matching_cols) {
    # If it exists, copy the values to the new "Cell_Type" field
    seuratObject@meta.data$Cell_Type <- seuratObject@meta.data[, col]
  }

  # Update the Seurat object in the list
  seurat_list[[i]] <- seuratObject
}

rm(seuratObject)
