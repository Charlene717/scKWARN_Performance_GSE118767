if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("forcats")) install.packages("forcats"); library(forcats)


plot_bar <- function(df, x, y, fill,
                     col.values = NULL,x_order = NULL, fill_order = NULL,
                     NameX = "X", NameY = "Y", NameFill = "Fill",
                     Set_TitleName="", aspect_ratio = 1) {

  # If x_order is not NULL, reorder x
  if (!is.null(x_order)) {
    df[[x]] <- factor(df[[x]], levels = x_order)
  }

  # If fill_order is not NULL, reorder fill
  if (!is.null(fill_order)) {
    df[[fill]] <- forcats::fct_relevel(df[[fill]], fill_order)
  }

  # Create the plot
  plot_Bar <- ggplot(df, aes(.data[[x]], .data[[y]], fill = .data[[fill]]))

  # Add the bar layer
  plot_Bar <- plot_Bar + geom_col(position = "dodge", alpha = 0.7)

  # Check the condition of col.values and fill_order
  if (!is.null(col.values) & !is.null(fill_order)) {
    plot_Bar <- plot_Bar + scale_fill_manual(values = col.values, breaks = fill_order)
  } else if (!is.null(col.values) & is.null(fill_order)) {
    plot_Bar <- plot_Bar + scale_fill_manual(values = col.values)
  } else if (is.null(col.values) & !is.null(fill_order)) {
    plot_Bar <- plot_Bar + scale_fill_discrete(breaks = fill_order)
  }

  # Set labels and theme
  plot_Bar <- plot_Bar +
    labs(title = Set_TitleName, x = NameX, y = NameY, fill = NameFill) +
    theme_bw() +
    theme(
      axis.line = element_line(linewidth = 0),  # set the size of the axis
      panel.border = element_rect(color = "black", linewidth = 1.2),  # bold border
      axis.title = element_text(size = 14, face = "bold"),  # enlarge the axis title
      axis.text = element_text(size = 14),  # enlarge the axis label
      axis.text.x = element_text(angle = 45, hjust = 1),  # adjust the angle of the x-axis label text
      legend.title = element_text(size = 14, face = "bold"),  # Enlarge and bold the legend title
      legend.text = element_text(size = 12),  # enlarge the legend text
      aspect.ratio = aspect_ratio
    ) +
    guides(fill = guide_legend(override.aes = list(size = 6)))  # Adjust the size of the points in the legend

  return(plot_Bar)
}



# ## Test function
# plot_bar(df = summary_df,
#          x = "Dataset",
#          y = "PCADepthCorr",
#          fill = "NorMeth",
#          col.values = col.values,
#          fill_order = Set_NorType,
#          NameX = "Dataset",
#          NameY = "Correlation",
#          NameFill = "Normalization",
#          Set_TitleName="Correlation Bar Plot")
