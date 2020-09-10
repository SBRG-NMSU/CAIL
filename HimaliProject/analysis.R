########### Prereqs ###########
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()
library(tidyverse)
library(xcms)

os <- Sys.info()
baseDir <- ifelse(os["sysname"] == "Windows", ifelse(os["nodename"] == "CAHF-RSCH-SK147",
    "C:/Users/ptrainor/Documents/GitHub/cail/","C:/Users/ptrainor/gdrive/CAIL/"), "~/gdrive/CAIL/")
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
iSTDs <- readxl::read_excel("iSTDs2.xlsx")

# Process profile annotation data:
profInfo <- as.data.frame(profInfo)
profInfo$profile_ID <- as.character(as.integer(profInfo$profile_ID))

# Import DDA data:
files1 <- list.files("C:/Users/ptrainor/Dropbox (NMSU Advanced)/Patrick/Himali/mzXML",
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
      l1[[j]]@intensity <- l1[[j]]@intensity / max(l1[[j]]@intensity) * 100
      l1[[j]] <- clean(removePeaks(l1[[j]], 2), all = TRUE)
      if(l1[[j]]@peaksCount == 0) l1[[j]] <- NULL
    }
    if(sum(sapply(l1, function(x) !is.null(x))) > 0){
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

save.image(file = "working_20200904.RData")

############ Modified Run MetFrag code ############
load("working_20200904.RData")
runMetFragMy <-function (config_file, MetFrag_dir, CL_name, config_dir = dirname(config_file)) {
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

save.image(file = "working_20200904b.RData")

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
presentISTD <- sampleAnno %>% select(ID, Label) %>% right_join(presentISTD, by = c("ID" = "file_ID"))

############ Some processing of profile data ############
# Set rownames on that data:
rownames(profs) <- sampleAnno$Label[match(rownames(profs), sampleAnno$ID)]

# Create data.frame with no blanks or "H20+IS samples":
aSProfs <- profs

# Which are 0:
which(aSProfs == 0, arr.ind = TRUE)

# Imputation Function
minNonZero <- function(x) min(x[x > 0]) / 2
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
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, numBins = 6), unit = "log2")
colEntropies <- apply(aSProfs, 2, entropyFun)
hist(colEntropies)
colEntropies <- as.data.frame(colEntropies)

# png(filename = "./Plots/EntropyFilter1.png", height = 3, width = 4, units = "in", res = 600)
ggplot(colEntropies, aes(x = colEntropies)) + geom_histogram(bins = 35, color = "black", fill = "lightblue") + 
  geom_vline(xintercept = 1, lty = 2, color = "darkred") + 
  labs(title = "Shannon entropy of peak intensities", x = "Entropy (Log2)", y = "Frequency") +
  theme_bw()
# dev.off()

############ Viz prior to normalization ############
# Heatmap #LOH
col1 <- colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256)

png(filename = "Plots/heatmap.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(aSProfs, colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

png(filename = "Plots/heatmap_wEntropy.png", height = 5, width = 7, units = "in", res = 300)
pheatmap::pheatmap(aSProfs[, colEntropies > 2], colsep = "", clustering_method = "ward.D2", color = col1,
                   show_colnames = FALSE)
dev.off()

# PCA
pca1 <- prcomp(aSProfs)
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
