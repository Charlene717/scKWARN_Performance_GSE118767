# 安裝並載入stringr 套件以使用 str_c 函數
if (!require("stringr")) {install.packages("stringr")}; library(stringr)

# 定義一個函數來處理時間的差異 (difftime 物件)
process_difftime <- function(time_diff) {
  time_value <- as.numeric(time_diff)
  time_unit <- attr(time_diff, "units")

  return(str_c(time_value, " ", time_unit))
}

# 將 difftime 物件轉換為適當的字串形式
Rec_Time_Spend.lt <- lapply(Rec_Time_Spend.lt, function(x) {
  if (inherits(x, "difftime")) {
    return(process_difftime(x))
  } else {
    return(as.character(x))
  }
})

# 將 POSIXct 物件轉換為適當的字串形式
Rec_Time_Point.lt <- lapply(Rec_Time_Point.lt, function(x) {
  if (inherits(x, "POSIXct")) {
    return(format(x, "%Y-%m-%d %H:%M:%S"))
  } else {
    return(as.character(x))
  }
})

# 將 Rec_Time_Spend.lt 和 Rec_Time_Point.lt 的內容分別轉換為文字
Rec_Time_Spend_Text <- str_c(names(Rec_Time_Spend.lt), Rec_Time_Spend.lt, sep = ": ", collapse = "\n")
Rec_Time_Point_Text <- str_c(names(Rec_Time_Point.lt), Rec_Time_Point.lt, sep = ": ", collapse = "\n")

# 組合所有文字，包含標題和空行
combined_text <- str_c("Rec_Time_Spend.lt:\n", Rec_Time_Spend_Text, "\n\nRec_Time_Point.lt:\n", Rec_Time_Point_Text)

# 將文字寫入到 txt 檔案
writeLines(combined_text, paste0(Name_ExportFolder,"/",Name_Export,"_Time_Record.txt"))
rm(combined_text)

#### Old version ####
# # 定義一個將列表資料寫入 txt 檔案的函數
# write_list_to_txt <- function(list_obj, file_name) {
#   # 將列表的內容轉換為文字
#   list_text <- str_c(names(list_obj), list_obj, sep = ": ", collapse = "\n")
#
#   # 將文字寫入到 txt 檔案
#   writeLines(list_text, file_name)
# }
#
# # 使用該函數將 Rec_Time_Spend.lt 和 Rec_Time_Point.lt 的內容寫入到 txt 檔案
# write_list_to_txt(Rec_Time_Spend.lt, paste0(Name_ExportFolder,"/",Name_Export,"_Rec_Time_Spend.txt"))
# write_list_to_txt(Rec_Time_Point.lt, paste0(Name_ExportFolder,"/",Name_Export,"_Rec_Time_Point.txt"))
