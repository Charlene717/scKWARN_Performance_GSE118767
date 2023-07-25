if (!require(ggplot2)) install.packages('ggplot2'); library(ggplot2)
if (!require(gridExtra)) install.packages('gridExtra'); library(gridExtra)


arrange_plots <- function(plot_list, ncol = 3, title_text = NULL, heights = NULL,
                          edgeOnlyXAxis = TRUE, edgeOnlyYAxis = TRUE,
                          edgeOnlyXLabel = TRUE, edgeOnlyYLabel = TRUE,
                          OneFigLegend = TRUE, removeTitles = FALSE,
                          x_tick_angle = 0, x_tick_size = 12,
                          y_tick_angle = 0, y_tick_size = 12,
                          removeX = FALSE, removeY = FALSE,auto_labels = TRUE,
                          legend_position = "bottom", legend_direction = "vertical",
                          simplifyXAxis = FALSE) {

  n <- length(plot_list)
  # Compute number of rows
  nrow <- ceiling(n / ncol)

  # Process X axis labels before modifying the plot_list
  if(simplifyXAxis){
    mapping <- list()
    for(i in 1:n){
      old_labels <- levels(plot_list[[i]]$data[[1]])  # Replace '1' with the correct index or name of your x variable
      new_labels <- 1:length(old_labels)
      mapping <- c(mapping, setNames(old_labels, new_labels))
      plot_list[[i]] <- plot_list[[i]] + scale_x_discrete(labels = new_labels)
    }

    mapping_df <- data.frame('Label' = names(mapping), 'Name' = unlist(mapping, use.names = FALSE))
    mapping_df$Label <- factor(mapping_df$Label, levels = rev(unique(mapping_df$Label))) # reverse order of labels
    mapping_plot <- ggplot(mapping_df, aes(x = Label, y = 1)) +
      geom_text(aes(label = paste(Label, ". ", Name)), hjust = 0) +
      theme_void() +
      coord_flip() +
      theme(plot.margin = margin(5, 40, 5, 5))  # Add some margin to the right

    plot_list <- c(plot_list, list(mapping_plot))
  }


  for(i in 1:n){
    # Get row and column index
    rowidx <- ceiling(i / ncol)
    colidx <- i %% ncol
    if(colidx == 0) colidx <- ncol

    # Modify axis.text (axis ticks)
    if(edgeOnlyXAxis && rowidx != nrow){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.text.x = element_blank())
    }
    if(edgeOnlyYAxis && colidx != 1){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.text.y = element_blank())
    }

    # Modify axis.title (axis labels)
    if(edgeOnlyXLabel && rowidx != nrow){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.title.x = element_blank())
    }
    if(edgeOnlyYLabel && colidx != 1){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.title.y = element_blank())
    }

    # Remove plot title if required
    if(removeTitles){
      plot_list[[i]] <- plot_list[[i]] + theme(plot.title = element_blank())
    }

    # Remove x-axis if required
    if(removeX){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
    }

    # Remove y-axis if required
    if(removeY){
      plot_list[[i]] <- plot_list[[i]] + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y =      element_blank())
    }

    # Adjust x and y axis tick text
    plot_list[[i]] <- plot_list[[i]] + theme(axis.text.x = element_text(size = x_tick_size, angle = x_tick_angle), axis.text.y = element_text(size = y_tick_size, angle = y_tick_angle))
  }

  # For single legend
  if(OneFigLegend){
    # Get legend from the first plot
    legend <- cowplot::get_legend(plot_list[[1]])
    # Remove legend from all plots
    for(i in 1:n){
      plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
    }
    # Create an empty plot to hold the legend
    legend_plot <- ggplot() + theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      annotation_custom(grob = legend)
    # Adjust the legend direction
    legend_plot <- legend_plot + guides(fill=guide_legend(direction = legend_direction))
    # Add it to the plot list
    plot_list <- c(plot_list, list(legend_plot))
  }

  # Combine all plots
  if(auto_labels) {
    labels <- "AUTO"
  } else {
    labels <- NULL
  }
  plot_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = ncol, align = "hv", axis = "tblr",
                                  labels = labels, hjust = -1, vjust = 1)
  # Add title if it's given
  if(!is.null(title_text)){
    title_plot <- ggplot() + theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0)) +
      ggtitle(title_text)
    plot_grid <- cowplot::plot_grid(title_plot, plot_grid, ncol = 1, rel_heights = c(0.1, 1))
  }

  # # For single legend
  # if(OneFigLegend){
  #   # Set the direction of legend before extracting it
  #   if (legend_direction == "horizontal") {
  #     plot_list[[1]] <- plot_list[[1]] + guides(fill=guide_legend(direction = legend_direction))
  #   }
  #
  #   # Get legend from the first plot
  #   legend <- cowplot::get_legend(plot_list[[1]])
  #
  #   # Remove legend from all plots
  #   for(i in 1:n){
  #     plot_list[[i]] <- plot_list[[i]] + theme(legend.position = "none")
  #   }
  #
  #   # Create an empty plot to hold the legend
  #   legend_plot <- ggplot() + theme_void() +
  #     theme(plot.margin = margin(0, 0, 0, 0)) +
  #     annotation_custom(grob = legend)
  #
  #   # Add it to the plot list
  #   plot_list <- c(plot_list, list(legend_plot))
  # }


  return(plot_grid)
}



# ## Test Function
# arrange_plots(plot_list_Bar)
# arrange_plots(plot_list_Bar[-1], edgeOnlyXAxis = TRUE, edgeOnlyYAxis = FALSE,
#               edgeOnlyXLabel = TRUE, edgeOnlyYLabel = FALSE, OneFigLegend = TRUE,
#               removeTitles = TRUE, legend_position = "bottom",x_tick_angle = 45,x_tick_size = 12,
#               title_text = "Test Title",removeX = TRUE,simplifyXAxis = TRUE) #,auto_labels = FALSE
#
# Set_SelectPlt <- c("ASWPCA_CellType","Purity","PCADepthCorr","ROGUE_CellType" , "ARI" )
# arrange_plots(plot_list_Bar[Set_SelectPlt], edgeOnlyXAxis = TRUE, edgeOnlyYAxis = FALSE,
#               edgeOnlyXLabel = TRUE, edgeOnlyYLabel = FALSE, OneFigLegend = TRUE,
#               removeTitles = TRUE, legend_position = "bottom",x_tick_angle = 45,x_tick_size = 12,
#               title_text = Name_FigTitle,removeX = TRUE,simplifyXAxis = TRUE,auto_labels = FALSE)
