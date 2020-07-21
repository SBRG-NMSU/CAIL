############ Prereqs ############
# Begin always run
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)
library(metfRag)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/ProducedWater"))
# End always run

############ Import MS/MS data ############
# Import: 
msmsData1 <- read.delim("MsmsIncluded_0_2020721724.txt", skip = 4, header = TRUE)

# Fix names:
names(msmsData1) <- gsub("\\.", "", names(msmsData1))

# Calculate Neutral mass:
msmsData1$neutralMass <- NA
msmsData1$neutralMass[msmsData1$Adducttype == "[M+H]+"] <- 
  msmsData1$AverageMz[msmsData1$Adducttype == "[M+H]+"] - 1.007276

msmsData1$neutralMass[msmsData1$Adducttype == "[M+Na]+"] <- 
  msmsData1$AverageMz[msmsData1$Adducttype == "[M+Na]+"] - 22.989218

msmsData1$neutralMass[msmsData1$Adducttype == "[M+NH4]+"] <- 
  msmsData1$AverageMz[msmsData1$Adducttype == "[M+NH4]+"] - 18.033823

# Functions:
getMSMS <- function(x){
  x2 <- str_split(x, " ", simplify = TRUE)[1,]
  x3 <- t(sapply(x2, FUN = function(y) str_split(y, "\\:", simplify = TRUE)))
  x3 <- apply(x3, 2, as.numeric)
  return(x3)
}
makeList <- function(x){
  rtMin <- x["AverageRtmin"]
  Mz <- x["AverageMz"]
  neutralMass <- x["neutralMass"]
  if(x["MSMSspectrum"] != ""){
    MSMS <- getMSMS(x["MSMSspectrum"])
  }else{
    MSMS <- ""
  }
  list(rtMin = as.numeric(rtMin), Mz = as.numeric(Mz), neutralMass = as.numeric(neutralMass), MSMS = MSMS)
}

# Call functions:
msmsData2 <- apply(msmsData1, 1, function(x) makeList(x))

############ Match MS1 ############
rtTol <- 20 / 60
mzTol <- 0.0005

test1 <- 9.253 > (msmsData1$AverageRtmin - rtTol) & 9.253 < (msmsData1$AverageRtmin + rtTol)
test2 <- 376.2602 > (msmsData1$AverageMz - mzTol) & 376.2602 < (msmsData1$AverageMz + mzTol)
which1 <- which(test1 & test2)

msmsData2[[which1]]$MSMS
msmsData2[[which1]]$neutralMass

############ MetFragR queries ############
sObj <- list()

sObj[["DatabaseSearchRelativeMassDeviation"]] <- 3.0
sObj[["FragmentPeakMatchAbsoluteMassDeviation"]] <- 0.001
sObj[["FragmentPeakMatchRelativeMassDeviation"]] <- 10
sObj[["PrecursorIonMode"]] <- 1
sObj[["IsPositiveIonMode"]] <- TRUE
sObj[["MetFragDatabaseType"]] <- "PubChem"
sObj[["NeutralPrecursorMass"]] <- 375.2526
sObj[["MetFragPreProcessingCandidateFilter"]] <- c("UnconnectedCompoundFilter","IsotopeFilter")
sObj[["MetFragPostProcessingCandidateFilter"]] <- c("InChIKeyFilter")
sObj[["PeakList"]] <- msmsData2[[which1]]$MSMS

cand1 <- run.metfrag(sObj)
