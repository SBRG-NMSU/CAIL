############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
library(tidyverse)
library(ggrepel)

oldPar <- par()
baseDir <- ifelse(Sys.info()["sysname"] == "Windows", "C:/Users/ptrainor/gdrive/CAIL/Energy/", "~/gdrive/CAIL/Energy/")
setwd(baseDir)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))