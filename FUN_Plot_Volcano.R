plot_volcano <- function(Marker.df,
                         x_field = "avg_log2FC",
                         y_field = "p_val",
                         gene_field = "Gene",
                         color = c(red = "#ef476f", gray = "gray", blue = "#0077b6"),
                         XThr = 1, YThr = 0.05,
                         ShowGeneNum = 5,
                         selected_genes = NULL,
                         Size_GeneLabel = 6,
                         x_limits = NULL,
                         y_limit = NULL){

  Xintercept = c(-XThr, XThr)
  Yintercept = -log10(YThr)

  library(ggplot2)
  library(cowplot)
  library(ggrepel)

  Marker.df2 <- Marker.df
  colnames(Marker.df2)[colnames(Marker.df2) == "Gene"] <- gene_field

  Marker.df2[[y_field]] <- Marker.df2[[y_field]] + 1.0e-300

  Marker.df2$color <- ifelse(Marker.df2[[y_field]] < YThr & abs(Marker.df2[[x_field]]) >= XThr, ifelse(Marker.df2[[x_field]] > XThr, 'red', 'blue'), 'gray')

  # Filter rows that pass thresholds
  Pos_df <- Marker.df2[Marker.df2[[x_field]] > XThr & Marker.df2[[y_field]] < YThr,]
  Neg_df <- Marker.df2[Marker.df2[[x_field]] < -XThr & Marker.df2[[y_field]] < YThr,]

  # Find the genes with highest and lowest x_field values
  Pos_List <- Pos_df[[gene_field]][order(-Pos_df[[x_field]])[1:ShowGeneNum]]
  Neg_List <- Neg_df[[gene_field]][order(Neg_df[[x_field]])[1:ShowGeneNum]]

  Marker.df2$genelabels <- factor(Marker.df2[[gene_field]], levels = c(Pos_List, Neg_List))

  # If selected_genes is provided, add these genes to the genelabels
  if(!is.null(selected_genes)){
    Marker.df2$genelabels <- factor(Marker.df2[[gene_field]], levels = unique(c(levels(Marker.df2$genelabels), selected_genes)))
    Marker.df2$selected <- ifelse(Marker.df2[[gene_field]] %in% selected_genes, "TRUE", "FALSE")
  } else {
    Marker.df2$selected <- "TRUE"
  }

  VolcanoPlot <- ggplot(Marker.df2, aes_string(x_field, paste0("-log10(", y_field, ")"), label = "genelabels")) +
    geom_point(aes(col = color), size = 3) +
    theme_bw() +
    labs(x = x_field, y = paste0("-log10(", y_field, ")")) +
    geom_hline(yintercept = Yintercept, lty = 8, col = "black", lwd = 0.8) +
    geom_vline(xintercept = Xintercept, lty = 8, col = "black", lwd = 0.8) +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          text = element_text(size = 15)) +
    geom_text_repel(aes(col = selected), na.rm = TRUE, size = Size_GeneLabel, box.padding = unit(0.45, "lines"), hjust = 1) +
    scale_color_manual(values = c("red" = "#ef476f", "gray" = "gray", "blue" = "#0077b6", "FALSE" = "#9ba0ab", "TRUE" = "#14213d")) +
    theme(aspect.ratio = 1)

  # Use coord_cartesian to adjust x and y axis limits
  if (!is.null(x_limits) || !is.null(y_limit)) {
    VolcanoPlot <- VolcanoPlot + coord_cartesian(xlim = c(-x_limits, x_limits), ylim = c(0, y_limit))
  }

  VolcanoPlot_2 <- VolcanoPlot + theme(panel.border = element_rect(fill = NA, color = "black", size = 2, linetype = "solid"))

  return(VolcanoPlot_2)
}

## Integrate Multiple VolcanoPlot
# source("FUN_Plot_Volcano.R")
plot_combined_volcanoes <- function(DEG_Input.list, Input_genes_list,
                                    Name_Dataset = "Dataset", Name_NorType = "LogNormalize",
                                    Name_ExportFolder = "ExportFolder", Name_Export = "Export",
                                    x_field = "logFC", y_field = "pvalue",
                                    gene_field = "primerid", XThr = 1, YThr = 0.05,
                                    Size_GeneLabel = 4, x_limits = NA, y_limit = NA,
                                    SubPlotTitle = "") {
  # Get the unique cluster IDs
  N_cluster <- names(DEG_Input.list)

  # Initialize a list to store the volcano plots
  volcano_plots <- list()

  # Loop through each cluster ID
  for (cluster_id in N_cluster) {
    tryCatch({
      # Get the result for the current cluster
      result <- DEG_Input.list[[cluster_id]]

      # Generate volcano plots
      if (is.na(x_limits) & is.na(y_limit)) {
        volcano_plot <- plot_volcano(result, x_field = x_field, y_field = y_field,
                                     gene_field = gene_field, XThr = XThr, YThr = YThr,
                                     Size_GeneLabel = Size_GeneLabel,
                                     selected_genes = Input_genes_list[[cluster_id]])
      } else {
        volcano_plot <- plot_volcano(result, x_field = x_field, y_field = y_field,
                                     gene_field = gene_field, XThr = XThr, YThr = YThr,
                                     Size_GeneLabel = Size_GeneLabel,
                                     selected_genes = Input_genes_list[[cluster_id]],
                                     x_limits = x_limits, y_limit = y_limit)
      }

      # Add title to the volcano plot
      volcano_plot <- volcano_plot + ggtitle(paste(SubPlotTitle, cluster_id))

      # Store the volcano plot
      volcano_plots[[cluster_id]] <- volcano_plot
    }, error = function(e) {
      warning(paste("Error in generating volcano plot for Cluster", cluster_id, ": ", conditionMessage(e)))
    })
  }

  # Combine volcano plots into a single plot
  if (!require("gridExtra")) install.packages("gridExtra")
  library(gridExtra)
  combined_volcano_plot <- do.call(grid.arrange, c(volcano_plots[order(as.numeric(names(volcano_plots)))], ncol = 3))

  library(grid)
  # Create custom title
  title_text <- textGrob(label = paste0(Name_Dataset, "_", Name_NorType), gp = gpar(fontsize = 16, fontface = "bold"))

  # Add title to the combined volcano plot
  combined_volcano_plot <- grid.arrange(title_text, combined_volcano_plot, ncol = 1, heights = c(0.05, 0.9))

  # Return the volcano plots
  return(volcano_plots)
}



# # Usage:
# source("Load_Marker_Gene.R")
# plot_volcano_Clt <- plot_combined_volcanoes(DEG_Cluster.list, selected_genes_list,
#                                             Name_Dataset = Name_Dataset, Name_NorType = Name_NorType,
#                                             Name_ExportFolder = Name_ExportFolder, Name_Export = Name_Export,
#                                             x_field = "logFC", y_field = "pvalue",
#                                             gene_field = "primerid", XThr = 1, YThr = 0.05,
#                                             Size_GeneLabel = 4,
#                                             x_limits = max(Result_Sum_Clt.df$logFC, na.rm = TRUE),
#                                             y_limit = max(-log10(Result_Sum_Clt.df$pvalue), na.rm = TRUE),
#                                             SubPlotTitle = "Cluster")
#
# plot_volcano_CT <- plot_combined_volcanoes(DEG_CellType.list, Marker_gene.list,
#                                            Name_Dataset = Name_Dataset, Name_NorType = Name_NorType,
#                                            Name_ExportFolder = Name_ExportFolder, Name_Export = Name_Export,
#                                            x_field = "logFC", y_field = "pvalue",
#                                            gene_field = "primerid", XThr = 1, YThr = 0.05,
#                                            Size_GeneLabel = 4,
#                                            x_limits = max(Result_Sum_CT.df$logFC, na.rm = TRUE),
#                                            y_limit = max(-log10(Result_Sum_CT.df$pvalue), na.rm = TRUE),
#                                            SubPlotTitle = "")
#

# # Access the DEG results for a specific cluster (e.g., cluster 1)
# DEG_CellType1 <- DEG_CellType.list[[1]]
#
#
# source("FUN_Plot_Volcano.R")
# plot_volcano(DEG_CellType.list[[1]], x_field = "logFC", y_field = "pvalue",gene_field = "primerid",
#                 XThr = 1, YThr = 0.05)
# plot_volcano(DEG_CellType.list[[1]], x_field = "logFC", y_field = "pvalue",gene_field = "primerid",
#                 XThr = 1, YThr = 0.05,selected_genes = c("CD14", "LYZ"))
# plot_volcano(DEG_CellType.list[[1]], x_field = "logFC", y_field = "pvalue",gene_field = "primerid",
#                 XThr = 1, YThr = 0.05,selected_genes = c("CD14", "LYZ"), x_limits = 4,
#                 y_limit = 300)
