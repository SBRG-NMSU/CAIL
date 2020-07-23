############ Prereqs ############
# Begin always run
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)
# library(metfRag)

# baseDir <- "C:/Users/ptrainor/gdrive/CAIL/"
baseDir <- "~/GitHub/cail/"
setwd(paste0(baseDir, "ProducedWater"))
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

############ Import Molecular Formula Enumeration ############


############ Match MS1 ############
rtTol <- 20 / 60
mzTol <- 3

# Add tolerances to MS1 data of MS/MS data:
msmsData1$rtLower <- msmsData1$AverageRtmin - rtTol
msmsData1$rtUpper <- msmsData1$AverageRtmin + rtTol
msmsData1$mzLower <- msmsData1$AverageMz - msmsData1$AverageMz * mzTol * 1e-6
msmsData1$mzUpper <- msmsData1$AverageMz + msmsData1$AverageMz * mzTol * 1e-6



############ MetFragR queries ############
for(i in 1:100){
  if(length(msmsData2[[i]]$MSMS) > 1){
    sObj <- list()
    cand <- NULL
    
    sObj[["DatabaseSearchRelativeMassDeviation"]] <- 3.0
    sObj[["FragmentPeakMatchAbsoluteMassDeviation"]] <- 0.1
    sObj[["FragmentPeakMatchRelativeMassDeviation"]] <- 10
    sObj[["PrecursorIonMode"]] <- 1
    sObj[["IsPositiveIonMode"]] <- TRUE
    sObj[["MetFragDatabaseType"]] <- "PubChem"
    sObj[["NeutralPrecursorMass"]] <- msmsData2[[i]]$neutralMass
    sObj[["MetFragPreProcessingCandidateFilter"]] <- c("UnconnectedCompoundFilter","IsotopeFilter")
    sObj[["MetFragPostProcessingCandidateFilter"]] <- c("InChIKeyFilter")
    sObj[["PeakList"]] <- msmsData2[[i]]$MSMS
    
    cand1 <- run.metfrag(sObj)
    
    # Check to see if any candidates were found; if so save:
    if(sum(cand1$Score) > 0){
      msmsData2[[i]]$candidates <- cand1
    }else{
      msmsData2[[i]]$candidates <- data.frame()
    }
    
  }else{
    msmsData2[[i]]$candidates <- data.frame()
  }
  print(i)
}

idk <- lapply(msmsData2, FUN = function(x) x$candidates)
