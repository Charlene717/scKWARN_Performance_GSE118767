if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)

create_box_plot <- function(df,
                            x_var = "", y_var = "",
                            fill_var = "",color_var = "",
                            Set_order = "No",
                            Set_Title = "", Name_XTitle = "", Name_YTitle = "",
                            x_text_size = 15, y_text_size = 15, legend_text_size = 12,
                            x_title_size = 17, y_title_size = 17,legend_title_size = 15,
                            color_values = NULL,
                            level_values = NULL,
                            legend_direction = "vertical") {
  # Preprocess data
  df[[y_var]] <- as.numeric(df[[y_var]])

  # Create initial plot
  if(Set_order == "sort"){
    plt.box <- ggplot(df,
                      aes_string(x = paste0("reorder(",x_var,",-",y_var,")"),
                                 y = y_var,
                                 fill = paste0("reorder(",fill_var,",-",y_var,")"),
                                 color = paste0("reorder(",color_var,",-",y_var,")")))
  } else if(Set_order == "revsort"){
    plt.box <- ggplot(df,
                      aes_string(x = paste0("reorder(",x_var,",",y_var,")"),
                                 y = y_var,
                                 fill = paste0("reorder(",fill_var,",",y_var,")"),
                                 color = paste0("reorder(",color_var,",",y_var,")")))
  } else if(Set_order == "user"){
    # Assign Order of Method
    if(!is.null(level_values)){
      df[[x_var]] <- factor(df[[x_var]], levels = level_values)
    }
    plt.box <- ggplot(df, aes_string(x = x_var, y = y_var, fill = fill_var, color = color_var))
  } else {
    plt.box <- ggplot(df, aes_string(x = x_var, y = y_var, fill = fill_var, color = color_var))
  }

  # Customize the plot
  plt.box <- plt.box  +
    geom_boxplot(alpha = 0, size = 1.2)

  # Check if custom colors are provided
  if(!is.null(color_values)){
    plt.box <- plt.box + scale_fill_manual(values = color_values) +
      scale_color_manual(values = color_values)
  }

  plt.box <- plt.box + theme_classic() +
    theme(
      panel.background = element_rect(fill = NA, colour = "black", size = 2),
      aspect.ratio = 1
    )

  ## Word Size
  plt.box <- plt.box + labs(x = Name_XTitle, y = Name_YTitle) +
    ggtitle(Set_Title) +
    labs(fill = Name_XTitle, color = Name_XTitle) +
    theme(
      title = element_text(size = 17),
      axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
      axis.title.x = element_text(size = x_title_size),# , face = "bold"),
      axis.text.y = element_text(size = y_text_size),
      axis.title.y = element_text(size = y_title_size),# , face = "bold"),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size)
    )

  # Add Mean to bar
  plt.box <- plt.box + stat_summary(fun = mean, geom = "text", aes(label = round(..y.., digits = 2)),
                                    vjust = -1.5, size = 6, color = "black")

  # Set the direction of the legend
  plt.box <- plt.box + guides(fill=guide_legend(direction = legend_direction), color=guide_legend(direction = legend_direction))

  return(plt.box)
}


# ## Test function
# NumCell <- ncol(Count.mtx)
# NumGene <- nrow(Count.mtx)
# # Results_Time.df <- read.csv("your_data.csv")  # replace "your_data.csv" with your actual file
# # Results_Time.df <- read.delim2(paste0(getwd(),"/Export/",ImportName,"_TimeRec.tsv"))
# plt.box <- create_box_plot(Results_Time.df,
#                            x_var = "Method", y_var = "TimeElapsed",
#                            fill_var = "Method",color_var = "Method",
#                            Set_order = "No",
#                            Set_Title = paste0(DatasetName," (",NumCell,"cells, ",NumGene,"genes)"),
#                            Name_XTitle = "Method",
#                            Name_YTitle = "Time spent(Sec)")
# plt.box
#
#
#
# # Set Color
# col.values <- c("counts" = "gray","scKWARN" = "blue", "RC" = "black", "scran" = "red", "SCNorm" = "#00CC00",
#                 "Seurat_LogNormalize" = "#159e5a","Seurat_RC" = "#2bed78","Seurat_CLR" = "#335e49",
#                 "sctransform" = "#B266FF", "PsiNorm" = "#FF8000")
# Set_NorType <- c("counts","scKWARN", "RC", "scran", "SCNorm",  # "SCNorm", # "Seurat_LogNormalize","Seurat_RC","Seurat_CLR",
#                  "sctransform", "PsiNorm")
#
# plt.box <- create_box_plot(Results_Time.df,
#                            x_var = "Method", y_var = "TimeElapsed",
#                            fill_var = "Method", color_var = "Method",
#                            Set_order = "user",
#                            Name_XTitle = "Method", Name_YTitle = "Time spent(Sec)",
#                            Set_Title = paste0(DatasetName," (",NumCell,"cells, ",NumGene,"genes)"),
#                            color_values = col.values,
#                            level_values = Set_NorType)
# plt.box


###################################################################################
create_box_plot_multiXTpye <- function(df,
                                       x_var = "", y_var = "",
                                       fill_var = "", color_var = "",
                                       Set_fill_order = NULL,
                                       Set_x_order = NULL,
                                       Set_Title = "", Name_XTitle = "", Name_YTitle = "", Name_LegendTitle = "",
                                       color_values = NULL, ylim_values = NULL,
                                       x_text_size = 15, y_text_size = 15, legend_text_size = 12,
                                       x_title_size = 17, y_title_size = 17,legend_title_size = 15,
                                       box_alpha = 0.2, jitter_alpha = 0.5,
                                       AspectRatio = 1, box_line = 1.2,rect_size = 1.3,
                                       show_mean = TRUE, legend_direction = "vertical") {

  # Preprocess data
  df[[y_var]] <- as.numeric(df[[y_var]])

  if(!is.null(Set_x_order)){
    df[[x_var]] <- factor(df[[x_var]], levels = Set_x_order)
  }

  # Check the order of fill_var
  if(!is.null(Set_fill_order)){
    df[[fill_var]] <- factor(df[[fill_var]], levels = Set_fill_order)
  }



  # Create initial plot
  plt.box <- ggplot(df, aes_string(x = x_var, y = y_var, fill = fill_var, color = color_var)) +
    geom_boxplot(position = position_dodge(0.9), na.rm = TRUE, alpha = box_alpha, size = box_line) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), alpha = jitter_alpha) +
    theme_classic() +
    theme_classic() +
    theme(
      panel.background = element_rect(fill = NA, colour = "black", size = rect_size),
      panel.grid.major.x = element_blank(),
      aspect.ratio = AspectRatio
    ) +
    labs(x = Name_XTitle, y = Name_YTitle) +
    ggtitle(Set_Title) +
    labs(fill = Name_LegendTitle, color = Name_LegendTitle) +
    theme(
      title = element_text(size = 17),
      axis.text.x = element_text(size = x_text_size, angle = 45, hjust = 1),
      axis.title.x = element_text(size = x_title_size),# , face = "bold"),
      axis.text.y = element_text(size = y_text_size),
      axis.title.y = element_text(size = y_title_size),# face = "bold"),
      legend.title = element_text(size = legend_title_size),
      legend.text = element_text(size = legend_text_size)
    )

  # Draw grey separators
  unique_x <- unique(df[[x_var]])
  for (i in 2:length(unique_x)) {
    plt.box <- plt.box + geom_vline(xintercept = i - 0.5, color = "grey", linetype = "dashed", size = 0.5)
  }

  # Check if custom colors are provided
  if(!is.null(color_values)){
    plt.box <- plt.box + scale_fill_manual(values = color_values) +
      scale_color_manual(values = color_values)

  }

  # Add Mean to bar if show_mean is TRUE
  if(show_mean){
    plt.box <- plt.box + stat_summary(
      fun = mean, na.rm = TRUE,
      geom = "text",
      aes(label = round(..y.., digits = 2), group = interaction(df[[x_var]], df[[fill_var]])),
      vjust = -2.5, size = 3, angle = 0,
      color = "black",
      position = position_dodge(0.9)
    )
  }

  # Set y axis limits if provided
  if(!is.null(ylim_values)){
    plt.box <- plt.box + ylim(ylim_values)
  }

  # Set the direction of the legend
  plt.box <- plt.box + guides(fill=guide_legend(direction = legend_direction), color=guide_legend(direction = legend_direction))


  return(plt.box)
}


# ## Test function
# Sum.df <- read.delim2(paste0(getwd(),"/Export/","MSINB_2023_06_05__scKW2000_Fea2000_SetCltNumTRUE_nPC50_FindVarFeadisp_V1_Sum_Test.txt"),
#                       na.strings = "NA")
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
#                                       Name_YTitle = "PCA Depth Correlation",
#                                       Name_LegendTitle = "Normalization",
#                                       Set_Title = "",
#                                       color_values = col.values,
#                                       ylim_values = c(0, 1), # set y limits
#                                       show_mean = TRUE) # show mean value on boxplot
# plt.box
#
# plt.box <- create_box_plot_multiXTpye(Sum.df,
#                                       x_var = "Dataset", y_var = "ROGUE_CellType",
#                                       fill_var = "NorMeth", color_var = "NorMeth",
#                                       Set_fill_order = Set_NorType,
#                                       Set_x_order = Set_Dataset,
#                                       Name_XTitle = "Dataset",
#                                       Name_YTitle = "ROGUE (CellType)",
#                                       Name_LegendTitle = "Normalization",
#                                       Set_Title = "",
#                                       color_values = col.values,
#                                       ylim_values = c(0, 1), # set y limits
#                                       show_mean = TRUE) # show mean value on boxplot
# plt.box

## To-Do List
# - [v] add separator line
# - [v] Order of Normalization
# - [v] x-axis ordering
# - [] Normalization affects the Dataset (x-axis) order
# - [v] Average position
# - [L] Aspect ratio adjustment
# - [v] ROGUE missing value handling
# - [v] Change the Function name
# - [v] Fig legend name correction

