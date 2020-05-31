############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/CAIL/WaterProject/AE_MS2_CmpID/", 
                  "~/gdrive/CAIL/WaterProject/AE_MS2_CmpID/")
setwd(baseDir)

library(tidyverse)

############ Import textfiles ############
# Annotation data from msp files:
load("C:/Users/ptrainor/gdrive/CAIL/WaterProject/AE_MS2_CmpID/mspInfo.RData")

# In Silico:
inSilico <- read.delim("inSilico_20200530.txt", header = TRUE, sep = "\t")
inSilico <- inSilico[, !names(inSilico) %in% c("File.path", "Title", "MS1.count", "X")]
colNames1 <- gsub("\\.1", "", colnames(inSilico)[grepl("\\.1", colnames(inSilico))])
colNames1b <- c("File.name", "MSMS.count", "PRECURSORMZ", "PRECURSORTYPE")

# Formula enumeration:
formEn <- read.delim("formulaEnumeration_20200530.txt", header = TRUE, sep = "\t")

############ Processing ############
# In Silico wide to long:
sList1 <- c("", 1:4)
inSilicoTemp2 <- list()
for(i in 1:length(sList1)){
  if(i > 1){
    colNames2 <- paste0(colNames1, ".", sList1[i])
  }else{
    colNames2 <- colNames1
    colNames2[colNames2 == "Structure.rank"] <- "Structure.rank.1"
  }
  
  inSilicoTemp <- inSilico[, (names(inSilico) %in% colNames1b | names(inSilico) %in% colNames2)]
  inSilicoTemp$rank <- i
  names(inSilicoTemp) <- gsub("\\.\\d", "", names(inSilicoTemp))
  inSilicoTemp2[[i]] <- inSilicoTemp
}
inSilico <- do.call("rbind", inSilicoTemp2)
rm(sList1, inSilicoTemp, inSilicoTemp2, colNames1, colNames1b, colNames2)

# Join with RT and precursor m/z:
df1$file <- gsub("\\.msp", "", df1$file)
inSilico <- df1 %>% left_join(inSilico, by = c("file" = "File.name"))
