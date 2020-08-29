########### Prereqs ###########
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()
library(tidyverse)
library(xcms)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/WaterProject"))

########### Import ###########
# Import DDA data:
files1 <- list.files("C:/Users/ptrainor/Dropbox (NMSU Advanced)/Patrick/Himali/mzXML",
                     full.names = TRUE, pattern = "MS2")
df1 <- readMSData(files1, mode = "onDisk")
df2 <- filterMsLevel(df1, 2)

precursorIntensity(df2)
