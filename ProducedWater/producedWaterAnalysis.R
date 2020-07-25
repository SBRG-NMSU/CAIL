############ Prereqs ############
# Begin always run
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)
library(doParallel)
# library(metfRag)

# baseDir <- "C:/Users/ptrainor/gdrive/CAIL/"
baseDir <- "~/GitHub/cail/"
setwd(paste0(baseDir, "ProducedWater"))
# End always run

############ Import enviMass profile data ############
load("enviMassOutput_20200721.RData")
profs0 <- profs
profs <- read.csv("profiled_peaks.csv")

profInfo2 <- profs[, !grepl("\\.sample\\.", names(profs))]
profs <- profs[, grepl("\\.sample\\.", names(profs))]

names(profs) <- gsub("\\.", "", gsub("X\\d\\.sample\\.", "", names(profs)))
profs <- t(profs)
colnames(profs) <- profInfo2$profile_ID

# Internal standards:
iSTDs <- readxl::read_excel("iSTDs.xlsx")

# Process profile annotation data:
profInfo <- as.data.frame(profInfo)
profInfo$profile_ID <- as.character(as.integer(profInfo$profile_ID))

############ Some processing of iSTD data ############
# Split ID from the enviMass output:
iSTDDF$iSTD_ID <- gsub("_none_none_none", "", iSTDDF$iSTD_ID)
temp1 <- str_split(iSTDDF$iSTD_ID, "_", simplify = TRUE)
iSTDDF$iSTD_ID <- temp1[,1]
iSTDDF$adduct <- temp1[,2]

# Join our db file:
iSTDs$ID <- as.character(iSTDs$ID)
iSTDDF <- iSTDDF %>% left_join(iSTDs, by = c("iSTD_ID" = "ID"))
iSTDDF$iso <- ifelse(abs(iSTDDF$`m/z` - iSTDDF$Da) < .050, "0", "")
mISTD <- iSTDDF %>% filter(adduct == "M+H" & iso == "0")
presentISTD <- mISTD %>% group_by(file_ID, Name) %>% summarize(n = n()) %>% spread(key = "Name", value = "n")

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

msmsData1$neutralMass[msmsData1$Adducttype == "[M+K]+"] <- 
  msmsData1$AverageMz[msmsData1$Adducttype == "[M+K]+"] - 38.963158

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
names(msmsData2) <- msmsData1$AlignmentID

############ Match MS1 ############
rtTol <- 5 / 60
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

# Add adduct information:
profInfo2$adduct <- gsub("\\*|\\*\\*", "", 
                         str_split(profInfo2$neutral_mass_info, " / ", simplify = TRUE)[,1])

# Aduct match table:
adductMatch <- data.frame(enviMassForm = c("M+H", "M+NH4", "M+Na", "M+K"), 
                          metFragAdduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+"))
profInfo2 <- profInfo2 %>% left_join(adductMatch, by = c("adduct"="enviMassForm"))

############ MetFrag queries ############
# Modified Run MetFrag code:
runMetFragMy <-function (config_file, MetFrag_dir, CL_name, config_dir = dirname(config_file)) {
  config_exists <- file.exists(config_file) && file.exists(config_dir)
  current_dir <- getwd()
  if (config_exists) {
    setwd(MetFrag_dir)
    log_dir <- gsub("config", "log", config_dir)
    if (!file.exists(log_dir)) {
      dir.create(log_dir)
    }
  }
  else {
    warning(paste("Configuration file ", config_file, 
                  " or directory not found, please try a new file"))
    stop()
  }
  MetFragCommand <- paste('java -jar "', gsub("/", "\\", paste0(MetFrag_dir, CL_name), fixed = TRUE), 
                          '" ', config_file, sep = "")
  MetFrag_out <- system(command = MetFragCommand, intern = TRUE, 
                        show.output.on.console = FALSE)
  log_file <- gsub("config", "log", config_file)
  write(MetFrag_out, log_file)
  setwd(current_dir)
}

# Base directory:
baseDir2 <- "C:/Users/ptrainor/Documents/GitHub/cail/ProducedWater"

profList <- profInfo2$profile_ID[1:10]

cl <- makeCluster(8)
registerDoParallel(cl)
foreach(i=1:10) %dopar%{
  profID <- profList[i]
  
  # Which MSMS from list object:
  whichMSMS <- profInfo2$msmsMatchU[profInfo2$profile_ID == profID]
  
  # Write temp .txt file with MS/MS data:
  write.table(msmsData2[[whichMSMS]]$MSMS, file = paste0(baseDir2, "/msmsPeaks/prof_", profID, ".txt"), 
              col.names = FALSE, row.names = FALSE, sep = "\t")
  
  # Write configuration file for MetFrag:
  profMZ <- profInfo2$profile_mean_mass[profInfo2$profile_ID == profID]
  profAdduct <- profInfo2$metFragAdduct[profInfo2$profile_ID == profID]
  ReSOLUTION::MetFragConfig(mass = profMZ, adduct_type = profAdduct, 
                            results_filename = paste0("prof_", profID),
                            peaklist_path = paste0(baseDir2, "/msmsPeaks/prof_", profID, ".txt"), 
                            base_dir = paste0(baseDir2, "/metFragOut"),
                            mzabs = 0.05, frag_ppm = 100, output = "CSV")
  
  # MetFrag call:
  runMetFragMy(paste0(baseDir2, "/metFragOut/config/prof_", profID, "_config.txt"), 
               MetFrag_dir = "C:/Program Files/metfrag/",
               CL_name = "MetFrag2.4.5-CL.jar")
}
stopCluster(cl)

baseDir <- "~/GitHub/cail/"
setwd(paste0(baseDir, "ProducedWater"))
