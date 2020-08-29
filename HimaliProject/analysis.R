########### Prereqs ###########
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()
library(tidyverse)
# library(xcms)

os <- Sys.info()
baseDir <- ifelse(os["sysname"] == "Windows", ifelse(os["nodename"] == "CAHF-RSCH-SK147","C:/Users/ptrainor/Documents/GitHub/cail/",
                                          "C:/Users/ptrainor/gdrive/CAIL/"), "~/gdrive/CAIL/")
setwd(paste0(baseDir, "HimaliProject/"))

########### Import and lite processing ###########
# Read in sample annotation:
sampleAnno <- read.csv("RData/sampleAnnotation_20200829.csv")
sampleAnno$ID <- as.character(sampleAnno$SampleNumber)

# Load profile data:
load("RData/enviMassOutput_20200826.RData")
rm(profs)
profs <- read.csv("RData/profiled_peaks_20200826.csv")
profInfo2 <- profs[, !grepl("X\\d", names(profs))]
profs <- profs[, grepl("X\\d", names(profs))]
profs <- t(profs)
colnames(profs) <- profInfo2$profile_ID
rownames(profs) <- gsub("X", "", str_split(rownames(profs), "\\.", simplify = TRUE)[,1])

# Process neutral mass / adduct data:
profInfo2$profile_ID <- as.character(profInfo2$profile_ID)
profInfo2$adduct <- gsub("\\*", "", str_split(profInfo2$neutral_mass_info, " / ", simplify = TRUE)[,1])

# Internal standards:
iSTDs <- readxl::read_excel("iSTDs.xlsx")

# Process profile annotation data:
profInfo <- as.data.frame(profInfo)
profInfo$profile_ID <- as.character(as.integer(profInfo$profile_ID))

# Import DDA data:
files1 <- list.files("C:/Users/ptrainor/Dropbox (NMSU Advanced)/Patrick/Himali/mzXML",
                     full.names = TRUE, pattern = "MS2")
df1 <- MSnbase::readMSData(files1, mode = "onDisk")
df2 <-  MSnbase::filterMsLevel(df1, 2)
df1Header <- MSnbase::header(df1)

########### Process MS/MS data ###########
df2Anno <- data.frame(mz = MSnbase::precursorMz(df2), precursorScan = MSnbase::precScanNum(df2))

rtTol <- 10 / 60
mzTol <- 2

# Add tolerances to MS1 data of MS/MS data:
msmsData1$rtLower <- msmsData1$AverageRtmin - rtTol
msmsData1$rtUpper <- msmsData1$AverageRtmin + rtTol
msmsData1$mzLower <- msmsData1$AverageMz - msmsData1$AverageMz * mzTol * 1e-6
msmsData1$mzUpper <- msmsData1$AverageMz + msmsData1$AverageMz * mzTol * 1e-6

# Join here:
profInfo2$msmsMatch <- profInfo2$msmsMatchU <- ""
for(i in 1:nrow(profInfo2)){
  match1 <- profInfo2$profile_mean_mass[i] > msmsData1$mzLower & profInfo2$profile_mean_mass[i] < msmsData1$mzUpper
  match2 <- profInfo2$profile_mean_RT_s[i] / 60 > msmsData1$rtLower & profInfo2$profile_mean_RT_s[i] / 60 < msmsData1$rtUpper
  match3 <- msmsData1$AlignmentID[which(match1 & match2)]
  
  profInfo2$msmsMatch[i] <- paste(match3, collapse = ";")
  
  if(length(match3) > 1){
    largeSN <- match3[which.max(msmsData1$SNaverage[msmsData1$AlignmentID %in% match3])]
    profInfo2$msmsMatchU[i] <- largeSN
  }
  else{
    profInfo2$msmsMatchU[i] <- profInfo2$msmsMatch[i]
  }
}
