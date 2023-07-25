plot_bubble <- function(df, x, y, size, fill,
                        col.values = NULL, x_order = NULL,fill_order = NULL,
                        Set_break = NULL, # Set_break = c(300, 1500, 3000),
                        NameX = "X", NameY = "Y",NameFill = "Fill",NameSize = "Size",
                        Set_TitleName="", aspect_ratio = 1,
                        fixed_size_value = 5,# fixed_size = FALSE,
                        Set_range = c(2, 10)) {

  # Ensure that column names exist in the data frame
  if (!(x %in% names(df) && y %in% names(df) && size %in% names(df) && fill %in% names(df))) {
    stop("One or more of the provided column names do not exist in the data frame.")
  }
  # # Ensure that column names are passed as strings
  # x <- deparse(substitute(x))
  # y <- deparse(substitute(y))
  # size <- deparse(substitute(size))
  # fill <- deparse(substitute(fill))
  #
  # # Check if the provided column names exist in the data frame
  # if (!all(c(x, y, size, fill) %in% names(df))) {
  #   stop("One or more of the provided column names do not exist in the data frame.")
  # }

    # Create the plot
    plot_Bubble <- ggplot(df, aes(!!sym(x), !!sym(y), fill = !!sym(fill)))

  # If x_order is not NULL, reorder x
  if (!is.null(x_order)) {
    plot_Bubble <- plot_Bubble + scale_x_discrete(limits = x_order)
  }

  if (is.null(Set_break)) {
    plot_Bubble <- plot_Bubble + geom_point(shape = 21, size = fixed_size_value, alpha = 0.7)
  } else {
    plot_Bubble <- plot_Bubble + geom_point(aes_string(size = size), shape = 21, alpha = 0.7) +
      scale_size_continuous(range = Set_range,breaks = Set_break)
  }

  plot_Bubble <- plot_Bubble + labs(title = Set_TitleName, color = "Norm")

  # Check the condition of col.values and fill_order
  if (!is.null(col.values) & !is.null(fill_order)) {
    plot_Bubble <- plot_Bubble + scale_fill_manual(values = col.values, breaks = fill_order)
  } else if (!is.null(col.values) & is.null(fill_order)) {
    plot_Bubble <- plot_Bubble + scale_fill_manual(values = col.values)
  } else if (is.null(col.values) & !is.null(fill_order)) {
    plot_Bubble <- plot_Bubble + scale_fill_discrete(breaks = fill_order)
  }


  plot_Bubble <- plot_Bubble +
    labs(x = NameX, y = NameY, size = NameSize, fill = NameFill) +
    theme_bw() +
    theme(
      axis.line = element_line(size = 0),  # set the size of the axis
      panel.border = element_rect(color = "black", size = 1.2),  # bold border
      axis.title = element_text(size = 14, face = "bold"),  # enlarge the axis title
      axis.text = element_text(size = 14),  # enlarge the axis label
      axis.text.x = element_text(angle = 45, hjust = 1),  # adjust the angle of the x-axis label text
      legend.title = element_text(size = 14, face = "bold"),  # Enlarge and bold the legend title
      legend.text = element_text(size = 12),  # enlarge the legend text
      aspect.ratio = aspect_ratio
    ) +
    guides(fill = guide_legend(override.aes = list(size = 6)))  # Adjust the size of the points in the legend

  return(plot_Bubble)
}


