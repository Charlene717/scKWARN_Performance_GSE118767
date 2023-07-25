

#######################
if(!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if(!require("Seurat")) install.packages("Seurat"); library(Seurat)
if(!require("ellipse")) install.packages("ellipse"); library(ellipse)

library(ggplot2)
library(dplyr)
library(purrr)



plot_clusters <- function(seuratObject, reduction = "umap", method = "none", cluster_label = "seurat_clusters") {
  if (!reduction %in% names(seuratObject@reductions)) {
    stop(paste0("The Seurat object does not contain ", reduction, ". Please add it before using this function."))
  }

  reduction_data <- as.data.frame(Embeddings(seuratObject, reduction))
  colnames(reduction_data) <- c(paste0(reduction, "_1"),
                                paste0(reduction, "_2"))

  reduction_data$cluster <- as.factor(seuratObject@meta.data[[cluster_label]])

  if (method == "none") {
    p <- ggplot(reduction_data, aes_string(x = names(reduction_data)[1], y = names(reduction_data)[2], color = "cluster")) +
      geom_point()

  } else if (method == "polygon") {
    convex_hulls <- reduction_data %>%
      split(.$cluster) %>%
      purrr::map(~.x[chull(.x[,1], .x[,2]), ])
    p <- ggplot(reduction_data, aes_string(x = names(reduction_data)[1], y = names(reduction_data)[2], color = "cluster")) +
      geom_point() +
      geom_polygon(data = do.call(rbind, convex_hulls), aes_string(fill = "cluster"), alpha = 0.2)

  } else if (method == "ellipse") {
    ellipse_data <- reduction_data %>%
      split(.$cluster) %>%
      purrr::map_df(function(df) {
        ellipse_points <- data.frame(ellipse::ellipse(cov(df[,1:2]), centre = colMeans(df[,1:2])))
        ellipse_points$cluster <- df$cluster[1]
        ellipse_points
      })
    names(ellipse_data)[1:2] <- names(reduction_data)[1:2]
    p <- ggplot(reduction_data, aes_string(x = names(reduction_data)[1], y = names(reduction_data)[2], color = "cluster")) +
      geom_point() +
      geom_path(data = ellipse_data, aes_string(x = names(reduction_data)[1], y = names(reduction_data)[2], color = "cluster"), size = 1)

  } else {
    stop(paste0("The method ", method , " is not recognized."))
  }
  p <- p + theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, size=1.5),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 17),
          legend.title = element_text(size = 17),
          legend.text = element_text(size = 15))


  p <- p +
    labs(# title = paste("Clustering of", reduction, "reduction"),
      # color = "Cluster",fill= "Cluster",
      x = paste0(reduction, "_1"),
      y = paste0(reduction, "_2")) +
    # scale_color_brewer(palette = "Set1") +
    theme(plot.title = element_text(hjust = 0.5, size = 20),
          legend.position = "top")+ theme(aspect.ratio=1)

  print(p)
}



# # Test function
#
# seuratObject <- seurat_Result_list[["mix.CELSeq51"]][["RC"]]
# seuratObject <- seurat_Result_list[["mix.CELSeq51"]][["scKWARN"]]
#
# plot_clusters(seuratObject, reduction = "pca", method = "none", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "umap", method = "none", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "tsne", method = "none", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "pca", method = "polygon", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "umap", method = "polygon", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "tsne", method = "polygon", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "pca", method = "ellipse", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "umap", method = "ellipse", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "tsne", method = "ellipse", cluster_label = "seurat_clusters")
# plot_clusters(seuratObject, reduction = "tsne", method = "ellipse", cluster_label = "Cell_Type")



################################################################################
# ####################### Old version #######################
# ##### polygon #####
#
# library(Seurat)
# # # 假设seuratObject是你的Seurat对象
# # seuratObject <- RunUMAP(seuratObject, dims = 1:10)
#
# # 现在我们来绘制UMAP，颜色按Cluster分类
# Seurat::DimPlot(seuratObject, group.by = "seurat_clusters")
#
#
# library(Seurat)
# library(ggplot2)
# library(dplyr)
# library(tidyverse)
#
# # 获取UMAP坐标和聚类结果
# umap_data <- as.data.frame(seuratObject@reductions$umap@cell.embeddings) %>%
#   rownames_to_column("cell") %>%
#   left_join(data.frame(cell = names(seuratObject@active.ident), cluster = seuratObject@active.ident), by = "cell")
#
# # 使用chull()函数计算每个簇的凸包
# convex_hulls <- umap_data %>%
#   split(.$cluster) %>%
#   map(~.x[chull(.x$UMAP_1, .x$UMAP_2), ])
#
# # 绘制UMAP图和凸包
# ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
#   geom_point() +
#   geom_polygon(data = do.call(rbind, convex_hulls), aes(fill = cluster), alpha = 0.2) +
#   theme_minimal()
#
#
#
#
# ##### ellipse #####
# # Install the ellipse package if you haven't done so
# if (!require(ellipse)) {
#   install.packages("ellipse")
# }
#
# library(ellipse)  # load the ellipse package
#
# library(ggplot2)
# library(dplyr)
# library(purrr)
#
#
# #### UMAP ####
# # Assuming seuratObject is your Seurat object
# seuratObject <- RunUMAP(seuratObject, dims = 1:10)
#
# umap_df <- as.data.frame(Embeddings(seuratObject, "umap"))
# colnames(umap_df) <- c("UMAP_1", "UMAP_2")
# umap_df$cluster <- as.factor(Idents(seuratObject))  # Convert cluster to factor
#
# # Calculate each cluster's ellipse centroid and generate ellipse data
# ellipse_data <- umap_df %>%
#   split(.$cluster) %>%
#   map_df(function(df) {
#     # Generate ellipse
#     ellipse_points <- data.frame(ellipse(cov(df[,1:2]), centre = colMeans(df[,1:2])))
#     # Add cluster column
#     ellipse_points$cluster <- df$cluster[1]
#     ellipse_points
#   })
#
# # Plot on UMAP with ellipse
# ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
#   geom_point() +
#   geom_path(data = ellipse_data, aes(x = UMAP_1, y = UMAP_2, color = cluster), size = 1) +
#   theme_minimal()
#
#
# ### PCA ####
#
# # Assuming seuratObject is your Seurat object
# seuratObject <- Runpca(seuratObject, dims = 1:10)
#
# pca_df <- as.data.frame(Embeddings(seuratObject, "pca"))
# colnames(pca_df) <- c("pca_1", "pca_2")
# pca_df$cluster <- as.factor(Idents(seuratObject))  # Convert cluster to factor
#
# # Calculate each cluster's ellipse centroid and generate ellipse data
# ellipse_data <- pca_df %>%
#   split(.$cluster) %>%
#   map_df(function(df) {
#     # Generate ellipse
#     ellipse_points <- data.frame(ellipse(cov(df[,1:2]), centre = colMeans(df[,1:2])))
#     # Add cluster column
#     ellipse_points$cluster <- df$cluster[1]
#     ellipse_points
#   })
#
# # Plot on pca with ellipse
# ggplot(pca_df, aes(x = pca_1, y = pca_2, color = cluster)) +
#   geom_point() +
#   geom_path(data = ellipse_data, aes(x = pca_1, y = pca_2, color = cluster), size = 1) +
#   theme_minimal()
#
#
# # #######################
# # ##### 很多輪廓線 #####
# # library(ggplot2)
# # library(ggExtra)
# # library(MASS)
# # library(purrr)
# # library(dplyr)
# #
# # # 假设seuratObject是你的Seurat对象
# # seuratObject <- RunUMAP(seuratObject, dims = 1:10)
# #
# # umap_df <- as.data.frame(Embeddings(seuratObject, "umap"))
# # colnames(umap_df) <- c("UMAP_1", "UMAP_2")
# # umap_df$cluster <- as.factor(Idents(seuratObject))  # Convert cluster to factor
# #
# # # 计算每个聚类的二维核密度估计，并取得等高线
# # contour_data <- umap_df %>%
# #   split(.$cluster) %>%
# #   map(function(df) {
# #     d <- kde2d(df$UMAP_1, df$UMAP_2, n = 100)
# #     contourLines(d$x, d$y, d$z)
# #   }) %>%
# #   map(function(lst) {
# #     do.call(rbind, lapply(lst, function(l) {
# #       data.frame(x = l$x, y = l$y, cluster = as.factor(l$level))  # Convert cluster to factor
# #     }))
# #   }) %>%
# #   bind_rows()
# #
# # # 在UMAP上绘制轮廓
# # ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
# #   geom_point() +
# #   geom_path(data = contour_data, aes(x = x, y = y, color = cluster), size = 1) +
# #   theme_minimal()
# #
