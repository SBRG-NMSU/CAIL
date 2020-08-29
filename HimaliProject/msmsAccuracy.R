########### Prereqs ###########
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()
library(tidyverse)

os <- Sys.info()
baseDir <- ifelse(os["sysname"] == "Windows", ifelse(os["nodename"] == "CAHF-RSCH-SK147","C:/Users/ptrainor/Documents/GitHub/cail/",
                                                     "C:/Users/ptrainor/gdrive/CAIL/"), "~/gdrive/CAIL/")
setwd(paste0(baseDir, "HimaliProject/"))

########### Load and process ###########
load("exampleMultSpec.RData")

df1 <- exampleMultSpec

df2 <- list()
for(i in 1:length(df1)){
  df2[[i]] <- data.frame(file = i, mz = df1[[i]]@mz, intensity = df1[[i]]@intensity)
  df2[[i]]$relIntensity <- df2[[i]]$intensity / sum(df2[[i]]$intensity)
}
df3 <- do.call("rbind", df2)
writexl::write_xlsx(df3, path = "combinedMSMS.xlsx")

