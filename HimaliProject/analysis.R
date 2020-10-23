########### Prereqs ###########
# Begin always run:
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()
library(tidyverse)
library(xcms)

os <- Sys.info()
baseDir <- ifelse(os["sysname"] == "Windows", ifelse(os["nodename"] == "CAHF-RSCH-SK147",
    "C:/Users/ptrainor/Documents/GitHub/cail/","C:/Users/ptrainor/gdrive/CAIL/"), "~/gdrive/CAIL/")
setwd(paste0(baseDir, "HimaliProject/"))
# End always run:

########### Import and lite processing ###########
# Read in sample annotation:
sampleAnno <- read.csv("RData/sampleAnnotation_20200829.csv")
sampleAnno$ID <- as.character(sampleAnno$SampleNumber)

# Load profile data:
load("RData/enviMassOutput_20201008.RData")
rm(profs)
profs <- read.csv("RData/profiled_peaks.csv")
profInfo2 <- profs[, !grepl("X\\d", names(profs))]
profs <- profs[, grepl("X\\d", names(profs))]
profs <- t(profs)
colnames(profs) <- profInfo2$profile_ID
rownames(profs) <- gsub("X", "", str_split(rownames(profs), "\\.", simplify = TRUE)[,1])

# Process neutral mass / adduct data:
profInfo2$profile_ID <- as.character(profInfo2$profile_ID)
profInfo2$adduct <- gsub("\\*", "", str_split(profInfo2$neutral_mass_info, " / ", simplify = TRUE)[,1])

# Internal standards:
iSTDs <- readxl::read_excel("iSTDs2.xlsx")

# Process profile annotation data:
profInfo <- as.data.frame(profInfo)
profInfo$profile_ID <- as.character(as.integer(profInfo$profile_ID))

# Import DDA data:
files1 <- list.files("C:/Users/ptrainor/Dropbox (NMSU Advanced)/Patrick/Himali2/mzXML",
                     full.names = TRUE, pattern = "MS2")
df1 <- readMSData(files1, mode = "onDisk")
df2 <-  filterMsLevel(df1, 2)
df1Header <- header(df1)
df1Header$sName <- rownames(df1Header)

########### Process MS/MS data ###########
df2Anno <- data.frame(precursorMz = MSnbase::precursorMz(df2), precursorScan = MSnbase::precScanNum(df2))
df2Anno$sName <- rownames(df2Anno)
df2Anno <- df2Anno %>% left_join(df1Header %>% select(sName, rt = retentionTime))

rtTol <- 10
mzTol <- 2

# Add tolerances to MS/MS data anno:
df2Anno$rtLower <- df2Anno$rt - rtTol
df2Anno$rtUpper <- df2Anno$rt + rtTol
df2Anno$mzLower <- df2Anno$precursorMz - df2Anno$precursorMz * mzTol * 1e-6
df2Anno$mzUpper <- df2Anno$precursorMz + df2Anno$precursorMz * mzTol * 1e-6

# Join here:
combFun <- function(x) consensusSpectrum(x, mzd = .1, minProp = .5, mzFun = mean, intensityFun = mean)
profInfo2$msmsMatch <- profInfo2$msmsMatchU <- ""
combSpectraList <- list()
for(i in 1:nrow(profInfo2)){
  match1 <- profInfo2$profile_mean_mass[i] > df2Anno$mzLower & profInfo2$profile_mean_mass[i] < df2Anno$mzUpper
  match2 <- profInfo2$profile_mean_RT_s[i] > df2Anno$rtLower & profInfo2$profile_mean_RT_s[i] < df2Anno$rtUpper
  match3 <- df2Anno$sName[which(match1 & match2)]
  
  if(length(match3) > 0){
    msms <- df2[rownames(df2@featureData) %in% match3]
    l1 <- list()
    for(j in 1:length(msms)){
      l1[[j]] <- msms[[j]]
      l1[[j]]<- clean(removePeaks(l1[[j]], 75), all = TRUE)
      if(l1[[j]]@peaksCount > 0){
        l1[[j]]@intensity <- l1[[j]]@intensity / max(l1[[j]]@intensity) * 100
        l1[[j]] <- clean(removePeaks(l1[[j]], 2), all = TRUE)
      }else{
        l1[[j]] <- NULL
      }
    }
    if(length(l1) > 0 && sum(sapply(l1, function(x) !is.null(x))) > 0){
      l1 <- rlist::list.clean(l1)
      mSpec1 <- MSpectra(l1)
      combMSMS <- combineSpectra(MSpectra(l1), method = combFun)
      combSpectraList[[profInfo2$profile_ID[i]]] <- combMSMS
      profInfo2$msmsMatch[i] <- paste(match3, collapse = ";")
    }
  }
  print(i)
}

# Add adduct information:
profInfo2$adduct <- gsub("\\*|\\*\\*", "", 
                         str_split(profInfo2$neutral_mass_info, " / ", simplify = TRUE)[,1])

# Aduct match table:
adductMatch <- data.frame(enviMassForm = c("M+H", "M+NH4", "M+Na", "M+K"), 
                          metFragAdduct = c("[M+H]+", "[M+NH4]+", "[M+Na]+", "[M+K]+"))
profInfo2 <- profInfo2 %>% left_join(adductMatch, by = c("adduct"="enviMassForm"))

save.image(file = "working_20201010.RData")

############ Modified Run MetFrag code ############
load("working_20201010.RData")
runMetFragMy <- function (config_file, MetFrag_dir, CL_name, config_dir = dirname(config_file)) {
  config_exists <- file.exists(config_file) && file.exists(config_dir)
  current_dir <- getwd()
  if (config_exists) {
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
                        show.output.on.console = TRUE)
  log_file <- gsub("config", "log", config_file)
  write(MetFrag_out, log_file)
  setwd(current_dir)
}

# Load EPA Comp-Tox file:
csvFile <- "C:/Users/ptrainor/Documents/GitHub/cail/ProducedWater/CompTox.csv"

# Run MetFrag with CompTox
profList <- profInfo2$profile_ID
for(i in 1:length(profList)){
  profID <- profList[i]
  
  # Which MSMS from list object:
  if(!is.null(combSpectraList[[profID]])){
    msms <- data.frame(mz = mz(combSpectraList[[profID]][[1]]), intensity = intensity(combSpectraList[[profID]][[1]]))
  }else{
    msms <- ""
  }
  
  # Write temp .txt file with MS/MS data:
  write.table(msms, file = paste0("msmsPeaks2/prof_", profID, ".txt"), 
              col.names = FALSE, row.names = FALSE, sep = "\t")
  
  # Write configuration file for MetFrag:
  profMZ <- profInfo2$profile_mean_mass[profInfo2$profile_ID == profID]
  profAdduct <- profInfo2$metFragAdduct[profInfo2$profile_ID == profID]
  ReSOLUTION::MetFragConfig(mass = profMZ, adduct_type = profAdduct, 
                            DB = "LocalCSV", localDB_path = csvFile, num_threads = 8,
                            results_filename = paste0("prof_", profID),
                            peaklist_path = paste0( "msmsPeaks2/prof_", profID, ".txt"), 
                            base_dir = paste0("metFragOut2"), ppm = 2,
                            mzabs = 0.1, frag_ppm = 333, output = "CSV")
  
  # MetFrag call:
  baseDir <- "C:/Users/ptrainor/Documents/GitHub/cail/HimaliProject"
  runMetFragMy(config_file = paste0(baseDir, "/metFragOut2/config/prof_", profID, "_config.txt"), 
               config_dir = paste0(baseDir, "/metFragOut2/config"), MetFrag_dir = "C:/Program Files/metfrag", 
               CL_name = "/MetFrag2.4.5-CL.jar")
  
  print(i)
}

save.image(file = "working_20201010b.RData")

############ Import MetFrag results ############
# Pickup here with previous save:
load("working_20201010b.RData")

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
  nrow1 <- nrow(res1)
  res1 <- res1 %>% filter(!grepl("Nb|\\[13C\\]|\\[14C\\]|D|Ag|Cu|Ba|Co|Zn|Ti|Hg|Pt|Si|Li|Fe|Re|Ir|Se|W|Ac|Al|Am|Sb|Ar|Be|Bi|B|Cd|Cs|Cf|Ce|Cr|Cm|CN|Es|Er|Eu|Gd|Ga|Ge|Au|Hf|Ho|In|Ir|Kr|La|Pb|Mg|Mn|Hg|Mo|Nd|Np|Ni|Os|Pd|Pu|Po|Pr|Pm|Ra|Rn|Re|Rh|Rb|Ru|Sm|Sc|Ta|Tc|Te|Tb|Tm|Sn|Ti|U|V|Xe|Yb|Y|Zn|Zr", res1$MolecularFormula))
  if(nrow(res1) != nrow1) print(i)
  
  if(nrow(res1) > 0){
    # Calculate mass error:
    res1$mzAbsErr <- res1$MonoisotopicMass - profInfo2$neutral_mass[profInfo2$profile_ID == profRes1]
    res1$mzRelErr <- abs(res1$mzAbsErr / res1$MonoisotopicMass * 1e6)
    
    # Add to all saved results:
    res1$prof <- profRes1
    res5[[i]] <- res1
    
    # If no MS/MS fragments then make new MS1 only res DF:
    if(is.null(combSpectraList[[profRes1]][[1]])){
      res2 <- res1[res1$mzRelErr == min(res1$mzRelErr),]
      res3 <- data.frame(prof = profRes1, formula = paste(unique(res2$MolecularFormula), collapse = "|"), 
                         id = "MS1_Only",
                         name = paste(paste(res2$CASRN_DTXSID, res2$MolecularFormula, res2$Name, sep = ":"), collapse = "|"),
                         InChI = "", InChIKey = "", SMILES = "")
      # Save top hit:
      res4[[i]] <- res3
      
    }else{
      # Get fragments:
      msmsRes1 <- data.frame(mz = mz(combSpectraList[[profRes1]][[1]]), intensity = intensity(combSpectraList[[profRes1]][[1]]))
      
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
        geom_segment(aes(x = `m/z`, xend = `m/z`, y = 0, yend = -value / 10), color = "#E41A1C", lwd = 1.2) + 
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

save(topHit, allHits, file = "metFragCand_20201010.RData")
save.image("working_20201010c.RData")

############ Some processing of iSTD data ############
load("working_20201010c.RData")

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
presentISTD <- sampleAnno %>% select(ID, Label) %>% right_join(presentISTD, by = c("ID" = "file_ID"))
presentISTD$Label[is.na(presentISTD$Label)] <- paste0("Pool", 1:5)

############ Some processing of profile data ############
# Set rownames on that data:
rownames(profs) <- sampleAnno$Label[match(rownames(profs), sampleAnno$ID)]
rownames(profs)[is.na(rownames(profs))] <- paste0("Pool", 1:5)
aSProfs <- profs

# Which are 0:
which(aSProfs == 0, arr.ind = TRUE)

# Remove those that are zero in all non-pool samples:
temp1 <- aSProfs[!grepl("Pool",rownames(aSProfs)),]
nonZeros <- apply(temp1, 2, function(x) sum(x > 0) > 0)
aSProfs <- aSProfs[, nonZeros]

# Imputation Function
minNonZero <- function(x) min(x[x > 0]) / 3
impFun <- function(x){
  if(sum(x == 0) > 0){
    x[x == 0] <- minNonZero(x)
  }
  return(x)
}
# Imputation for missing values:
aSProfs <- apply(aSProfs, 2, impFun)

# Convert to log-scale:
aSProfs <- log10(aSProfs)

# Entropy filter:
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, numBins = 8), unit = "log2")
colEntropies <- apply(aSProfs[!grepl("Pool",rownames(aSProfs)),], 2, entropyFun)
hist(colEntropies)
colEntropies <- as.data.frame(colEntropies)

png(filename = "./Plots/Entropy.png", height = 4, width = 5.5, units = "in", res = 300)
ggplot(colEntropies, aes(x = colEntropies)) + geom_histogram(bins = 20, color = "black", fill = "lightblue") + 
  geom_vline(xintercept = 1.25, lty = 2, color = "darkred") + 
  labs(title = "Shannon entropy of peak intensities", x = "Entropy (Log2)", y = "Frequency") +
  theme_bw()
dev.off()

############ Viz prior to normalization ############
# Heatmap
col1 <- rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))

png(filename = "Plots/heatmap.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(aSProfs, colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

png(filename = "Plots/heatmap_wEntropy.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(aSProfs[, colEntropies > 1.25], colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

# Heatmap on centered data:
png(filename = "Plots/heatmapCentered.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(scale(aSProfs, scale = FALSE), colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

png(filename = "Plots/heatmapCentered_wEntropy.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(scale(aSProfs[, colEntropies > 1.25], scale = FALSE), colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

# Heatmap on centered and scaled data:
png(filename = "Plots/heatmapCenteredScaled.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(scale(aSProfs, scale = TRUE), colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

png(filename = "Plots/heatmapCenteredScaled_wEntropy.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(scale(aSProfs[, colEntropies > 1.25], scale = TRUE), colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

# PCA prior to norm
pca1 <- prcomp(aSProfs[!grepl("Pool",rownames(aSProfs)), colEntropies > 1.25])
pca1DF <- as.data.frame(pca1$x[, 1:4])
pca1DF$sampleID <- rownames(pca1DF)

png(filename = paste0("./Plots/PCA_prior2Norm_", gsub("-", "", Sys.Date()), ".png"),
    height = 5.5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(pca1DF, aes(x = PC1, y = PC2, label = sampleID)) + geom_point() + 
  ggrepel::geom_text_repel(size = 2.75) + theme_bw() + labs(title = "PCA prior to normalization") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(filename = paste0("./Plots/PCA_prior2Norm_23_", gsub("-", "", Sys.Date()), ".png"),
    height = 5.5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(pca1DF, aes(x = PC2, y = PC3, label = sampleID)) + geom_point() + 
  ggrepel::geom_text_repel(size = 2.75) + theme_bw() + labs(title = "PCA prior to normalization") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# PCA prior to norm scaled
pca1 <- prcomp(aSProfs[!grepl("Pool",rownames(aSProfs)), colEntropies > 1.25], scale = TRUE)
pca1DF <- as.data.frame(pca1$x[, 1:4])
pca1DF$sampleID <- rownames(pca1DF)

png(filename = paste0("./Plots/PCA_prior2Norm_Scaled_", gsub("-", "", Sys.Date()), ".png"),
    height = 5.5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(pca1DF, aes(x = PC1, y = PC2, label = sampleID)) + geom_point() + 
  ggrepel::geom_text_repel(size = 2.75) + theme_bw() + labs(title = "PCA prior to normalization (scaled)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

png(filename = paste0("./Plots/PCA_prior2Norm_23_Scaled_", gsub("-", "", Sys.Date()), ".png"),
    height = 5.5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(pca1DF, aes(x = PC2, y = PC3, label = sampleID)) + geom_point() + 
  ggrepel::geom_text_repel(size = 2.75) + theme_bw() + labs(title = "PCA prior to normalization (scaled)") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

########### Intensity distributions ###########
# Profile data wide to long:
profs2 <- profs
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName)
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
# Filter out missing:
profs2b <- profs2 %>% filter(intensity > 0)
p1 <- ggplot(profs2b, aes(x = fileName, y = log10(intensity))) + geom_boxplot() + 
  geom_hline(yintercept = median(log10(profs2b$intensity)), color = "darkblue", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

# How many profiles found?
profs3 <- as.data.frame(profs)
profs3$fileName <- rownames(profs3)
profs3$fileName <- factor(profs3$fileName)
profs3 <- profs3 %>% gather(key = "profID", value = "intensity", -fileName)
profs3 <- profs3 %>% mutate(isFound = ifelse(intensity > 0, 1, 0))
profs3 <- profs3 %>% group_by(fileName) %>% summarize(profiles = sum(isFound))

p2 <- ggplot(profs3, aes(x = fileName, y = profiles)) + geom_point() + 
  geom_hline(yintercept = median(profs3$profiles), color = "darkred", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

png(filename = paste0("./Plots/Intens_",gsub("-", "", Sys.Date()), ".png"),
    height = 5, width = 7, units = "in", res = 600)
p1
dev.off()

png(filename = paste0("./Plots/PeakCount_",gsub("-", "", Sys.Date()), ".png"),
    height = 5, width = 7, units = "in", res = 600)
p2
dev.off()

png(filename = paste0("./Plots/IntensPeakCount_",gsub("-", "", Sys.Date()), ".png"),
    height = 10, width = 7, units = "in", res = 600)
gridExtra::grid.arrange(p1, p2, nrow = 2)
dev.off()

########### Median deviation normalization ###########
medByFile <- profs2b %>% group_by(fileName) %>% summarize(medInt = median(intensity))
medNorm <- medByFile %>% mutate(grandMed = median(medInt), medRatio = medInt / grandMed, multFactorDist = 1 / medRatio)

########### IS-based normalization ###########
mISTD2 <- mISTD %>% group_by(file_ID, Name) %>% select(file_ID, Name, Intensity) %>% as.data.frame()
mISTD2 <- mISTD2 %>% left_join(presentISTD %>% select(ID, Label), by = c("file_ID" = "ID"))

# Remove one outlier:
mISTD2 <- mISTD2 %>% filter(!(Name == "Atrazine-D5" & file_ID == 16) & 
                              !(Name == "Caffeine-13C3" & file_ID == 1))
iSTDCV <- mISTD2 %>% group_by(Name) %>% summarize(cv = sd(Intensity) / mean(Intensity))

# Filter for the ones that work:
mISTD3 <- mISTD2 %>% filter(Name %in% iSTDCV$Name[iSTDCV$cv < .20])
mISTD3 <- mISTD3 %>% group_by(Name) %>% 
  mutate(medInt = median(Intensity, na.rm = TRUE), medRatio = Intensity / medInt)
mISTD4 <- mISTD3 %>% group_by(Label) %>% summarize(medMedRatio = median(medRatio), 
                                                   multFactorISTD = 1 / medMedRatio) 

p4 <- ggplot(mISTD3, aes(x = Name, y = log10(Intensity), color = Label, group = Label)) + geom_point() + 
  geom_line() + theme_bw() + labs(color = "Type") + theme(axis.text.x = element_text(angle = 90))

png(filename = paste0("./Plots/iSTDIntens_",gsub("-", "", Sys.Date()), ".png"),
    height = 4, width = 8, units = "in", res = 600)
p4
dev.off()

# Join two normalization df's:
medNorm <- medNorm %>% left_join(mISTD4, by = c("fileName" = "Label"))
medNorm$multFactor <- .5 * (medNorm$multFactorISTD + medNorm$multFactorDist)
medNorm$multFactor <- .5 * (1 + medNorm$multFactor)

# Scale the data:
sProfs <- profs
for(i in 1:nrow(sProfs)){
  sProfs[i, ] <- sProfs[i, ] * medNorm$multFactor[match(rownames(sProfs)[i], medNorm$fileName)]
}

# Make a new boxplot
profs2 <- sProfs
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName)
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
# Filter out missing:
profs2b <- profs2 %>% filter(intensity > 0)
p3 <- ggplot(profs2b, aes(x = fileName, y = log10(intensity))) + geom_boxplot() + 
  geom_hline(yintercept = median(log10(profs2b$intensity)), color = "darkblue", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

png(filename = paste0("./Plots/IntensNorm_",gsub("-", "", Sys.Date()), ".png"),
    height = 5, width = 7, units = "in", res = 600)
p3
dev.off()

rm(colEntropies, iSTDCV, iSTDDF, iSTDs, medNorm, mISTD, mISTD2, mISTD3, mISTD4, p1, p2, p3, p4, p5,
   pca1, pca1DF, presentISTD, temp1, entropyFun, medByFile, profs2b, profs3)

############ Get classyfire classification for top hits ############
# Import MetFrag results:
hitsList <- list()
for(i in 1:nrow(topHit)){
  if(topHit$InChI[i] != ""){
    hit <- classyfireR::get_classification(topHit$InChIKey[i])
    if(!is.null(hit)){
      hitClassDF <- hit@classification
      hitClassDF$Level <- factor(hitClassDF$Level, levels = hitClassDF$Level)
      hitClassDF <- hitClassDF %>% select(-CHEMONT) %>% spread(key = Level, value = Classification)
      hitClassDF$prof <- topHit$prof[i]
      hitsList[[i]] <- hitClassDF
      cat("Successful: ", i, "\n")
    }
    else{
      cat("Not successful: ", i, "\n")
    }
  }
}
hitList <- do.call("bind_rows", hitsList)
hitList <- hitList %>% select(prof, kingdom, superclass, class, subclass, level5 = `level 5`, 
                              level6 = `level 6`, level7 = `level 7`, level8 = `level 8`)
topHit <- topHit %>% left_join(hitList)

save.image("working_20201010d.RData")

############ Export ############
load("working_20201010d.RData")

sProfsT <- as.data.frame(t(sProfs))
sProfsT$id <- rownames(sProfsT)
exportDF <- profInfo2 %>% left_join(topHit, by = c("profile_ID" = "prof")) %>% 
  filter(id == "MS1 + MS2") %>% 
  select(id = profile_ID, mass = profile_mean_mass, rtSec = profile_mean_RT_s, 
         rtMin = profile_mean_RT_min, adduct = metFragAdduct, formula, name, superclass, 
         class, subclass, subsubclass = level5, subsubsubclass = level6) %>%
  left_join(sProfsT)
writexl::write_xlsx(exportDF, path = "export.xlsx")
