FUN_TimeRec <- function(InputData = Count.mtx,
                        Method_Name = "Nor_Seurat_Log",
                        Set_Times = 20,
                        TimePoint.lt = list(), TimeSpend.lt = list(),
                        Time.df = data.frame(Method = character(), Times = integer(), TimeElapsed = numeric()),
                        func = Seurat::NormalizeData, ...){
  for (i in 1:Set_Times) {
    TimePoint.lt[[paste0(Method_Name,"_S",i)]] <- Sys.time()

    func(InputData,...)

    TimePoint.lt[[paste0(Method_Name,"_E",i)]] <- Sys.time()
    TimeSpend.lt[[paste0(Method_Name,i)]] <- as.numeric(TimePoint.lt[[paste0(Method_Name,"_E",i)]] - TimePoint.lt[[paste0(Method_Name,"_S",i)]],
                                                                   units = "secs")

    # Add results to data.frame
    Time.df <- Time.df %>% add_row(Method = Method_Name , Times = i, TimeElapsed = as.numeric(TimeSpend.lt[[paste0(Method_Name,i)]]))

  }

  OutPut <- list(TimePoint.lt = TimePoint.lt,
                 TimeSpend.lt = TimeSpend.lt,
                 Time.df = Time.df
  )
  return(OutPut)
}

# ##### Example ######
# set.seed(12345)
# G <- 2000; n <- 600 # G: number of genes, n: number of cells
# NB_cell <- function(j) rnbinom(G, size = 0.1, mu = rgamma(G, shape = 2, rate = 2))
# countsimdata <- sapply(1:n, NB_cell)
# colnames(countsimdata) <- paste("cell", 1:n, sep = "_")
# rownames(countsimdata) <- paste("gene", 1:G, sep = "_")
# ## Test FUN_TimeRec
# Test.lt <- FUN_TimeRec(countsimdata)
# Rec_Time_Point.lt <- Test.lt[["TimePoint.lt"]]; Rec_Time_Spend.lt <- Test.lt[["TimeSpend.lt"]]; Results_Time.df <- Test.lt[["Time.df"]]
