## Function of cell number Plot

if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("gridExtra")) install.packages("gridExtra"); library(gridExtra)

# Define a function to create a bar plot for each Seurat object
create_plots_CellCount_list <- function(seurat_list, reverse_order = FALSE,
                                        add_sum_cells = FALSE, YLimCellNum = NULL, Set_Cell = "Cell_Type") {
  # Reverse the list order if required
  if (reverse_order) {
    seurat_list <- rev(seurat_list)
  }

  # Create a list of plots
  plots_CellCount_list <- lapply(names(seurat_list), function(name) {
    df <- seurat_list[[name]]@misc[["CountCell"]][[Set_Cell]]

    # Add sum of cells to the title if required
    title <- name
    if (add_sum_cells) {
      sum_cells <- sum(df$count)
      title <- paste0(name, " (Sum Cells: ", sum_cells, ")")
      seurat_list[[name]]@misc[["CountCell"]][["Sum"]] <- sum_cells
    }

    # Create a histogram
    plot <- ggplot(df, aes_string(x = Set_Cell, y = "count", fill = Set_Cell)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_fill_discrete() +
      labs(title = title, x = NULL, y = "Count") +
      theme_bw() +
      theme(panel.grid = element_line(color = "#ededed"),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5),
            aspect.ratio = 1,
            axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 13),
            axis.title.y = element_text(size = 12)) +
      theme(panel.border = element_rect(linetype = "solid", colour = "black", size = 1.1)) +
      geom_text(aes_string(label = "count", color = Set_Cell), vjust = -0.5, size = 5)

    # Check the order of YLimCellNum
    if(!is.null(YLimCellNum)){
      plot <- plot + coord_cartesian(ylim = c(0, YLimCellNum))
      # plot <- plot + ylim(0, YLimCellNum)
    }

    return(plot)
  })

  return(plots_CellCount_list)
}

