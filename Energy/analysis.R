############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
library(tidyverse)

oldPar <- par()
baseDir <- ifelse(Sys.info()["sysname"] == "Windows", "C:/Users/ptrainor/gdrive/CAIL/Energy/", "~/gdrive/CAIL/Energy/")
setwd(baseDir)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process ############
df1 <- readxl::read_excel("data_20200914.xlsx")
table(df1$`sample ID`)
