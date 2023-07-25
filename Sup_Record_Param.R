# 將參數轉換成一個字串
param_str <- paste("Set_logPlusOne:", Set_logPlusOne, "\n",
                   "Set_logDepth:", Set_logDepth, "\n",
                   "Set_numGeneforEst:", Set_numGeneforEst, "\n",
                   "Set_nfeatures:", Set_nfeatures, "\n",
                   "Num_PCA:", Num_PCA, "\n",
                   "Set_NorType:", paste(Set_NorType, collapse = ", "), "\n",
                   "Set_sctransform_Res:", Set_sctransform_Res, "\n",
                   "Set_Metrics:", paste(Set_Metrics, collapse = ", "), "\n",
                   "Set_Dataset:", paste(Set_Dataset, collapse = ", "), "\n",
                   "Name_CP:", Name_CP, "\n",
                   "Name_Sup:", Name_Sup, "\n",
                   "Name_Note:", Name_Note, "\n",
                   "Name_FigTitle:", Name_FigTitle, "\n",
                   "Name_Export:", Name_Export)

# 將字串寫入文件
writeLines(param_str, con = paste0(Name_ExportFolder,"/",Name_Export,"_parameters.txt"))



