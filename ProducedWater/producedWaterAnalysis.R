############ Prereqs ############
# Begin always run
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)

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
  if(is.null(dim(x3))) x3 <- matrix(x3, ncol = 2, byrow = TRUE)
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

# Add adduct information:
profInfo2$adduct <- gsub("\\*|\\*\\*", "", 
                         str_split(profInfo2$neutral_mass_info, " / ", simplify = TRUE)[,1])

# Aduct match table:
adductMatch <- data.frame(enviMassForm = c("M+H", "M+NH4", "M+Na", "M+K"), 
                          metFragAdduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+"))
profInfo2 <- profInfo2 %>% left_join(adductMatch, by = c("adduct"="enviMassForm"))

############ Modified Run MetFrag code ############
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

############ EPA ChemTox MetFrag queries ############
compToxFile <- "C:/Users/ptrainor/Dropbox (NMSU Advanced)/Patrick/DSSTOX_MS_Ready_Chemical_Structures/DSSToxMS-Ready.xlsx"
compTox <- readxl::read_excel(compToxFile)

# Match to example .csv format:
compToxMapper <- data.frame(currentNames = c("Preferred_Name", "CAS-RN", "DSSTox_Substance_ID", "SMILES", 
                            "SMILES (MS-Ready)", "Formula", "Monoisotopic Mass (MS-Ready)",
                            "InChIString (MS-Ready)", "InChIKey (MS-Ready)", "DSSTox_Compound_ID (MS-Ready)"),
           newNames = c("Name", "CASRN_DTXSID", "Identifier", "SMILES_DTXSID", 
                        "SMILES", "MolecularFormula", "MonoisotopicMass",
                        "InChI", "InChIKey", "DTXCID_MSready"))
compTox <- compTox %>% select(compToxMapper$currentNames)
names(compTox) <- compToxMapper$newNames[match(names(compTox), compToxMapper$currentNames)]

# Remove incomplete cases:
isNC <- rep(TRUE, nrow(compTox))
for(i in 1:nrow(compTox)){
  isNC[i] <- !complete.cases(compTox[i, ])
}
compTox <- compTox[!isNC,]

# Break InChIKey into pieces:
InChIKeyTemp <- str_split(compTox$InChIKey, "-", simplify = TRUE)
compTox$InChIKey1 <- InChIKeyTemp[,1]
compTox$InChIKey2 <- InChIKeyTemp[,2]
compTox$InChIKey3 <- InChIKeyTemp[,3]

# Remove "Aux" info from InChIString:
compTox$InChI <- str_split(compTox$InChI, "\n", simplify = TRUE)[,1]

csvFile <- "C:/Users/ptrainor/Documents/GitHub/cail/ProducedWater/CompTox.csv"
#csvFile <- "C:/Users/ptrainor/Documents/GitHub/cail/ProducedWater/kegg_2017-07-23.csv"

write.csv(compTox, csvFile, na = "", row.names = FALSE)

# Base directory:
baseDir2 <- "C:/Users/ptrainor/Documents/GitHub/cail/ProducedWater"

# Stop here and save:
save.image("working_20200817.RData")

# Run MetFrag with CompTox
profList <- profInfo2$profile_ID
for(i in 1:length(profList)){
  profID <- profList[i]
  
  # Which MSMS from list object:
  whichMSMS <- profInfo2$msmsMatchU[profInfo2$profile_ID == profID]
  
  # Write temp .txt file with MS/MS data:
  write.table(msmsData2[[whichMSMS]]$MSMS, file = paste0(baseDir2, "/msmsPeaks2/prof_", profID, ".txt"), 
              col.names = FALSE, row.names = FALSE, sep = "\t")
  
  # Write configuration file for MetFrag:
  profMZ <- profInfo2$profile_mean_mass[profInfo2$profile_ID == profID]
  profAdduct <- profInfo2$metFragAdduct[profInfo2$profile_ID == profID]
  ReSOLUTION::MetFragConfig(mass = profMZ, adduct_type = profAdduct, 
                            DB = "LocalCSV", localDB_path = csvFile,
                            results_filename = paste0("prof_", profID),
                            peaklist_path = paste0(baseDir2, "/msmsPeaks2/prof_", profID, ".txt"), 
                            base_dir = paste0(baseDir2, "/metFragOut2"), ppm = 2,
                            mzabs = 0.1, frag_ppm = 333, output = "CSV")
  
  # MetFrag call:
  runMetFragMy(paste0(baseDir2, "/metFragOut2/config/prof_", profID, "_config.txt"), 
               MetFrag_dir = "C:/Program Files/metfrag/", CL_name = "MetFrag2.4.5-CL.jar")
  
  print(i)
}

############ Pubmed online MetFrag queries ############
# profList <- profInfo2$profile_ID
# for(i in 1:length(profList)){
#   profID <- profList[i]
#   
#   # Which MSMS from list object:
#   whichMSMS <- profInfo2$msmsMatchU[profInfo2$profile_ID == profID]
#   
#   # Write temp .txt file with MS/MS data:
#   write.table(msmsData2[[whichMSMS]]$MSMS, file = paste0(baseDir2, "/msmsPeaks/prof_", profID, ".txt"), 
#               col.names = FALSE, row.names = FALSE, sep = "\t")
#   
#   # Write configuration file for MetFrag:
#   profMZ <- profInfo2$profile_mean_mass[profInfo2$profile_ID == profID]
#   profAdduct <- profInfo2$metFragAdduct[profInfo2$profile_ID == profID]
#   ReSOLUTION::MetFragConfig(mass = profMZ, adduct_type = profAdduct, 
#                             results_filename = paste0("prof_", profID),
#                             peaklist_path = paste0(baseDir2, "/msmsPeaks/prof_", profID, ".txt"), 
#                             base_dir = paste0(baseDir2, "/metFragOut"), ppm = 2,
#                             mzabs = 0.05, frag_ppm = 100, output = "CSV")
#   
#   # MetFrag call:
#   runMetFragMy(paste0(baseDir2, "/metFragOut/config/prof_", profID, "_config.txt"), 
#                MetFrag_dir = "C:/Program Files/metfrag/", CL_name = "MetFrag2.4.5-CL.jar")
#   
#   print(i)
# }

############ Import MetFrag results ############
# Pickup here with previous save:
load("working_20200817.RData")

# Get all the file names from the directory:
fNames1 <- list.files('metFragOut2/results/')

# List for storing everything:
res5 <- list()
# List for storing top hit:
res4 <- list()

for(i in 1:length(fNames1)){
  fNameRes1 <- fNames1[i]
  profRes1 <- str_split(fNameRes1, "_|\\.", simplify = TRUE)[,2]
  res1 <- read.csv(paste0("metFragOut2/results/", fNameRes1))
  
  # Filter for super not sensible elements:
  res1 <- res1 %>% filter(!grepl("Nb|\\[13C\\]|\\[14C\\]|D|Ag|Cu|Ba|Co|Zn|Ti|Hg|Pt|Si|Li|Fe|Re|Ir|Se|W|Ac|Al|Am|Sb|Ar|Be|Bi|B|Cd|Cs|Cf|Ce|Cr|Cm|CN|Es|Er|Eu|Gd|Ga|Ge|Au|Hf|Ho|In|Ir|Kr|La|Pb|Mg|Mn|Hg|Mo|Nd|Np|Ni|Os|Pd|Pu|Po|Pr|Pm|Ra|Rn|Re|Rh|Rb|Ru|Sm|Sc|Ta|Tc|Te|Tb|Tm|Sn|Ti|U|V|Xe|Yb|Y|Zn|Zr", res1$MolecularFormula))
  
  if(nrow(res1) > 0){
    # Calculate mass error:
    res1$mzAbsErr <- res1$MonoisotopicMass - profInfo2$neutral_mass[profInfo2$profile_ID == profRes1]
    res1$mzRelErr <- abs(res1$mzAbsErr / res1$MonoisotopicMass * 1e6)
    
    # Add to all saved results:
    res1$prof <- profRes1
    res5[[i]] <- res1
    
    # Get fragments:
    msmsRes1 <- msmsData2[[profInfo2$msmsMatchU[profInfo2$profile_ID == profRes1]]]$MSMS
    
    # If no MS/MS fragments then make new MS1 only res DF:
    if(is.null(msmsRes1) || msmsRes1 == ""){
      res2 <- res1[res1$mzRelErr == min(res1$mzRelErr),]
      res3 <- data.frame(prof = profRes1, formula = paste(unique(res2$MolecularFormula), collapse = "|"), 
                         id = "MS1_Only",
                         name = paste(paste(res2$CASRN_DTXSID, res2$MolecularFormula, res2$Name, sep = ":"), collapse = "|"),
                         InChI = "", InChIKey = "", SMILES = "")
      # Save top hit:
      res4[[i]] <- res3
      
    }else{
      # Make sure these are ordered by score:
      res1 <- res1 %>% arrange(desc(Score))
      
      # Save top result:
      res2 <- res1[1,]
      res3 <- data.frame(prof = profRes1, formula = res2$MolecularFormula, id = "MS1 + MS2",
                         name = paste(res2$CASRN_DTXSID, res2$Name, sep = ":"),
                         InChI = res2$InChI, InChIKey = paste0(res2$InChIKey1, "-", res2$InChIKey2, "-", res2$InChIKey3),
                         SMILES = res2$SMILES)
      # Save top hit:
      res4[[i]] <- res3
      
      ## Plot the top scoring:
      # Fragment match:
      msmsRes1 <- as.data.frame(msmsRes1)
      names(msmsRes1) <- c("m/z", "Intensity")
      msmsRes1$Matched <- "No"
      msmsRes1$value <- NA
      
      # Get MetFrag match results:
      if(!is.na(res1$ExplPeaks[1])){
        expPeaksRes1 <- str_split(str_split(res1$ExplPeaks[1], ";", simplify = TRUE), "_", simplify = TRUE)
        expPeaksRes1 <- as.matrix(expPeaksRes1)
        expPeaksRes1 <- as.data.frame(expPeaksRes1)
        for(j in 1:ncol(expPeaksRes1)) expPeaksRes1[,j] <- as.numeric(expPeaksRes1[,j])
        names(expPeaksRes1)[1] <- "m/z"
        
        # Match fragments and MetFrag attributed:
        for(j in 1:nrow(msmsRes1)){
          which1 <- which(msmsRes1$`m/z`[j] < expPeaksRes1$`m/z` + .0001 & 
                            msmsRes1$`m/z`[j] > expPeaksRes1$`m/z` - .0001)
          if(length(which1) > 0){
            msmsRes1$Matched[j] <- "Yes"
            msmsRes1$value[j] <- expPeaksRes1$V2[which1]
          }
        }
      }
      
      # Make a spectral match plot:
      tempInfo <- profInfo2[profInfo2$profile_ID == profRes1,]
      png(file = paste0("metFragOut2/specMatchPlots/", "profMatch_", tempInfo$profile_ID, ".png"),
          height = 5, width = 7, units = "in", res = 300)
      p1 <- ggplot(msmsRes1, aes(x = `m/z`, xend = `m/z`, y = 0, yend = Intensity, color = Matched)) + 
        geom_segment(lwd = 1.2) + theme_bw() + scale_color_brewer(palette = "Set1") +
        labs(y = "Intensity", title = paste0("Profile: ", tempInfo$profile_ID, ". precursor m/z: ", 
                                             round(tempInfo$profile_mean_mass,4), ", RT: ", tempInfo$profile_mean_RT_min, " min"),
             subtitle = paste0("Candidate. ", res2$MolecularFormula, ", CAS: ", res2$CASRN_DTXSID,
                               ", Name: ", res2$Name))
      show(p1)
      dev.off()
      
      png(file = paste0("metFragOut2/specMatchPlots/", "profScore_", tempInfo$profile_ID, ".png"),
          height = 7, width = 7, units = "in", res = 300)
      p2 <- ggplot(msmsRes1, aes(x = `m/z`, xend = `m/z`, y = 0, yend = Intensity)) + 
        geom_segment(lwd = 1.2, color = "#377EB8") + 
        geom_segment(aes(x = `m/z`, xend = `m/z`, y = 0, yend = -value), color = "#E41A1C", lwd = 1.2) + 
        geom_hline(yintercept = 0, lwd = .1) + theme_bw() + 
        labs(y = "Intensity (-Fragment Score)", title = paste0("Profile: ", tempInfo$profile_ID, ". precursor m/z: ", 
                                                               round(tempInfo$profile_mean_mass,4), ", RT: ", tempInfo$profile_mean_RT_min, " min"),
             subtitle = paste0("Candidate. ", res2$MolecularFormula, ", CAS: ", res2$CASRN_DTXSID,
                               ", Name: ", res2$Name)) 
      show(p2)
      dev.off()
    }
  }
  print(i)
}
res4 <- do.call("rbind", res4)
res5 <- do.call("rbind", res5)

topHit <- res4
allHits <- res5

save(topHit, allHits, file = "metFragCand_20200818.RData")

############ Join with profile data ############
profsL <- as.data.frame(t(profs))
profsL$prof <- as.integer(rownames(profsL))
res4$prof <- as.integer(res4$prof)
profsL <- profInfo2 %>% select(prof = profile_ID, rt = profile_mean_RT_s, rtMin = profile_mean_RT_min, mz = profile_mean_mass,
                               neutral_mass, adduct = metFragAdduct, ISTD, homologue) %>% full_join(res4) %>% full_join(profsL)
profsL <- profsL %>% filter(ISTD == "-")
profsL <- profsL %>% select(-ISTD)
writexl::write_xlsx(profsL, path = "ProducedWater_20200821.xlsx")
