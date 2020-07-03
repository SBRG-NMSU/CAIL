############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)
library(ggrepel)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/WaterProject"))

############ Functions ############
# Flag all missing:
misCheckFun <- function(x){
  if(sum(x == 0) == length(x)){
    logical1 <- TRUE
  }else{
    logical1 <- FALSE
  }
  return(logical1)
}

misCheckFun2 <- function(x){
  if(sum(x == 0) > .10 * length(x)){
    logical1 <- TRUE
  }else{
    logical1 <- FALSE
  }
  return(logical1)
}

# Min non zero half function:
minNonZero <- function(x) min(x[x > 0]) / 2

# Imputation Function
impFun <- function(x){
  if(sum(x == 0) > 0){
    x[x == 0] <- minNonZero(x)
  }
  return(x)
}

############ Imports ############
# Read in sample annotation:
sampleAnno <- read.csv("SampleAnno_20200514.csv")
sampleAnno$ID <- as.character(sampleAnno$ID)
sampleAnno$Name <- factor(sampleAnno$Name, levels = sampleAnno$Name[sampleAnno$myOrder])

# Load profile data:
load("enviMassOutput_20200509.RData")

# Internal standards:
iSTDs <- readxl::read_excel("iSTDs.xlsx")

# Process profile annotation data:
profInfo <- as.data.frame(profInfo)
profInfo$profile_ID <- as.character(as.integer(profInfo$profile_ID))

# Import ID information:
fenIn <- readxl::read_excel("AE_MS2_CmpID/AE_CmpdID_20200620.xlsx", sheet = "FormulaEn_inSilico")
fenIn$feature <- as.character(as.integer(fenIn$feature))
# Unique:
fenIn2 <- fenIn %>% select(feature, precursorMZ, RT) %>% unique()

# Import MS/MS data:
msmsData <- readxl::read_excel("AE_MS2_CmpID/AE_CmpdID_20200620.xlsx", sheet = "peakData")
msmsData$feature <- as.character(msmsData$feature)

############ Some processing of MS/MS data ############
msmsDataTemp1 <- str_split(msmsData$MS2Spectrum, " ")
procFun0 <- function(x){
  temp1 <- str_split(x, ":", simplify = TRUE)
  temp1 <- as.data.frame(temp1)
  names(temp1) <- c("mz", "intensity")
  temp1 <- temp1 %>% filter(intensity > 100)
  return(temp1)
}
msmsData2 <- lapply(msmsDataTemp1, procFun0)
names(msmsData2) <- msmsData$feature
rm(msmsDataTemp1, msmsDataTemp2)

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
presentISTD <- sampleAnno %>% select(ID, Name, Type) %>% right_join(presentISTD, by = c("ID" = "file_ID"))

############ Some processing of profile data ############
# Set rownames on that data:
rownames(profs) <- sampleAnno$Name[match(rownames(profs), sampleAnno$ID)]

# Create data.frame with no blanks or "H20+IS samples":
actSamps <- sampleAnno$Name[sampleAnno$Type == "sample" & sampleAnno$profiled == TRUE]
aSProfs <- profs[rownames(profs) %in% actSamps,]

# Which are 0:
which(aSProfs == 0, arr.ind = TRUE)

# Which are not missing in the samples:
aSProfs <- aSProfs[, !apply(aSProfs, 2, misCheckFun)]

# Imputation for missing values:
aSProfs <- apply(aSProfs, 2, impFun)

# Convert to log-scale:
aSProfs <- log2(aSProfs)
aSProfs <- scale(aSProfs, center = TRUE, scale = FALSE)

# Entropy filter:
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, 5), unit = "log2")
colEntropies <- apply(aSProfs, 2, entropyFun)
hist(colEntropies)
aSProfs <- aSProfs[, colEntropies > 1]
colEntropies <- as.data.frame(colEntropies)

# png(filename = "./Plots/EntropyFilter1.png", height = 3, width = 4, units = "in", res = 600)
ggplot(colEntropies, aes(x = colEntropies)) + geom_histogram(bins = 35, color = "black", fill = "lightblue") + 
  geom_vline(xintercept = 1, lty = 2, color = "darkred") + 
  labs(title = "Shannon entropy of peak intensities", x = "Entropy (Log2)", y = "Frequency") +
  theme_bw()
# dev.off()

############ Viz prior to normalization ############
# Heatmap
heatmap(aSProfs)

# PCA
pca1 <- prcomp(aSProfs)
pca1DF <- as.data.frame(pca1$x[, 1:2])
pca1DF$sampleID <- rownames(pca1DF)
pca1DF$sampleID <- gsub("90percent_", "", gsub("90per_", "", pca1DF$sampleID))
pca1DF <- pca1DF %>% left_join(sampleAnno %>% select(Name, Group = tag3), by = c("sampleID" = "Name"))

# png(filename = "./Plots/PCA_prior2Norm.png", height = 5.5, width = 8, units = "in", res = 600)
set.seed(3)
ggplot(pca1DF, aes(x = PC1, y = PC2, label = sampleID, color = Group)) + geom_point() + 
  geom_text_repel(size = 2.75) + theme_bw() + labs(title = "PCA prior to normalization") +
  theme(plot.title = element_text(hjust = 0.5))
# dev.off()

########### Intensity distributions ###########
# Profile data wide to long:
profs2 <- profs
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName, levels = levels(sampleAnno$Name))
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
# Filter out missing:
profs2b <- profs2 %>% filter(intensity > 0)
p1 <- ggplot(profs2b, aes(x = fileName, y = log10(intensity))) + geom_boxplot() + 
  geom_hline(yintercept = median(log10(profs2b$intensity)), color = "darkblue", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

# How many profiles found?
profs3 <- as.data.frame(profs)
profs3$fileName <- rownames(profs3)
profs3$fileName <- factor(profs3$fileName, levels = levels(sampleAnno$Name))
profs3 <- profs3 %>% gather(key = "profID", value = "intensity", -fileName)
profs3 <- profs3 %>% mutate(isFound = ifelse(intensity > 0, 1, 0))
profs3 <- profs3 %>% group_by(fileName) %>% summarize(profiles = sum(isFound))

p2 <- ggplot(profs3, aes(x = fileName, y = profiles)) + geom_point() + 
  geom_hline(yintercept = median(profs3$profiles), color = "darkred", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

# png(filename = paste0("./Plots/Intens_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 5, width = 7, units = "in", res = 600)
p1
# dev.off()

# png(filename = paste0("./Plots/PeakCount_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 5, width = 7, units = "in", res = 600)
p2
# dev.off()

# png(filename = paste0("./Plots/IntensPeakCount_",gsub("-", "", Sys.Date()), ".png"),
#     height = 10, width = 7, units = "in", res = 600)
gridExtra::grid.arrange(p1, p2, nrow = 2)
# dev.off()

########### Median deviation normalization ###########
medByFile <- profs2b %>% group_by(fileName) %>% summarize(medInt = median(intensity))
medNorm <- medByFile %>% mutate(grandMed = median(medInt), medRatio = medInt / grandMed, multFactorDist = 1 / medRatio)

########### IS-based normalization ###########
mISTD2 <- mISTD %>% group_by(file_ID, Name) %>% select(file_ID, Name, Intensity) %>% as.data.frame()
mISTD2 <- mISTD2 %>% left_join(sampleAnno %>% select(ID, fileName = Name, grp = tag3), by = c("file_ID" = "ID"))
iSTDCV <- mISTD2 %>% group_by(Name, grp) %>% summarize(cv = sd(Intensity) / mean(Intensity))

# Filter for the four that work (and remove fake QC samples):
mISTD3 <- mISTD2 %>% filter(Name %in% c(c("Acetaminophen-D4","Atenolol-D7", "Atrazine-D5", "Caffeine-13C3")))
mISTD3 <- mISTD3 %>% filter(grp != "FALSE")
mISTD3 <- mISTD3 %>% group_by(Name) %>% 
  mutate(medInt = median(Intensity, na.rm = TRUE), medRatio = Intensity / medInt)
mISTD4 <- mISTD3 %>% group_by(fileName) %>% summarize(medMedRatio = median(medRatio), multFactorISTD = 1 / medMedRatio) 

p4 <- ggplot(mISTD3, aes(x = Name, y = log10(Intensity), color = grp, group = fileName)) + geom_point() + 
  geom_line() + theme_bw() + labs(color = "Type")

# png(filename = paste0("./Plots/iSTDIntens_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 5, width = 7, units = "in", res = 600)
p4
# dev.off()

# Plot of scaling:
# png(filename = paste0("./Plots/iSTDsf_",gsub("-", "", Sys.Date()), ".png"),
#      height = 5, width = 6, units = "in", res = 600)
ggplot(mISTD4, aes(x = fileName, y = log10(medMedRatio))) + geom_point() + 
  geom_hline(yintercept = 0, lty = 2, color = "darkred") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Current Scale") 
# dev.off()

# From manual integration:
mISTD4$medMedRatio[mISTD4$fileName == "RO_3_2hr"] <- 10^(-.25)
mISTD4$multFactorISTD[mISTD4$fileName == "RO_3_2hr"] <- 1 / mISTD4$medMedRatio[mISTD4$fileName == "RO_3_2hr"]

# Join two normalization df's:
medNorm <- medNorm %>% left_join(mISTD4)

# iSTD inverse factor:
ggplot(mISTD4, aes(x = fileName, y = log10(medMedRatio))) + geom_point() + 
  geom_hline(yintercept = 0, lty = 2, color = "darkred") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y = "Current Scale") 

# Scale the data:
medNorm$multFactor <- 10^((log10(medNorm$multFactorISTD) + 3*log10(medNorm$multFactorDist)) / 4)
sProfs <- profs
for(i in 1:nrow(sProfs)){
  sProfs[i, ] <- sProfs[i, ] * medNorm$multFactor[match(rownames(sProfs)[i], medNorm$fileName)]
}

# Scale the internal standards:
mISTD3 <- mISTD3 %>% left_join(medNorm, by = "fileName") %>% mutate(Intensity2 = Intensity * multFactor)
p5 <- ggplot(mISTD3, aes(x = Name, y = log10(Intensity2), color = grp, group = fileName)) + geom_point() + 
  geom_line() + theme_bw() + labs(color = "Type")
# png(filename = "./Plots/iSTDIntens2.png", height = 9, width = 7, units = "in", res = 600)
gridExtra::grid.arrange(p4, p5, nrow = 2)
# dev.off()

# Make a new boxplot
profs2 <- sProfs
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName, levels = levels(sampleAnno$Name))
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
# Filter out missing:
profs2b <- profs2 %>% filter(intensity > 0)
p3 <- ggplot(profs2b, aes(x = fileName, y = log10(intensity))) + geom_boxplot() + 
  geom_hline(yintercept = median(log10(profs2b$intensity)), color = "darkblue", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "")

gridExtra::grid.arrange(p1, p3, nrow = 2)

rm(colEntropies, iSTDCV, iSTDDF, iSTDs, medNorm, mISTD, mISTD2, mISTD3, mISTD4, p1, p2, p3, p4, p5,
   pca1, pca1DF, presentISTD, temp1, actSamps, entropyFun, misCheckFun, medByFile, profs2b, profs3)

save.image("working_20200702.RData")

########### Algal time course ###########
load("working_20200702.RData")

# First get data for the profiles for AE:
profs3 <- profs2 %>% filter(grepl("AE-", fileName))
profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- ifelse(grepl("AE-3|AE-4", profs3$fileName), "AE-3/4", "AE-1")

# Processing for time-course:
# 1. Remove profiles not observed in >40% of AE-1, AE-3, and AE-4 samples
# 2. Remove profiles with average peak area less than 5×10^6 in both AE-1 and AE-3/4

# Filter out super missing:
profs3 %>% select(profID) %>% unique() %>% nrow() # 48,255 profiles
profs3 <- profs3 %>% group_by(profID) %>% mutate(countNonZero = sum(intensity > 0), countTotal = n(),
                                                 removeIt = ifelse(countNonZero / countTotal <= .6, TRUE, FALSE)) %>%
  ungroup() %>% filter(!removeIt) %>% select(-(countNonZero:removeIt))
profs3 %>% select(profID) %>% unique() %>% nrow() # 1,933 profiles

# Filter out too low abundance: 
profs3 <- profs3 %>% group_by(profID, whichSP) %>% mutate(meanIntByGroup = mean(intensity), 
     lowAbundanceByGroup = meanIntByGroup < 5*10^6, groupN = n()) %>% ungroup() %>% group_by(profID) %>% 
     mutate(countLowAbundanceByGroup = sum(lowAbundanceByGroup), tooLow = countLowAbundanceByGroup == groupN) %>%
  group_by(profID) %>% mutate(removeThis = sum(tooLow) > 0) %>% ungroup() %>% filter(!removeThis) %>%
  select(-meanIntByGroup, -lowAbundanceByGroup, -groupN, -countLowAbundanceByGroup, -tooLow, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 1,460 profiles

# Imputation
profs3 <- profs3 %>% group_by(profID) %>% mutate(mNZ = minNonZero(intensity))
profs3$intensity <- ifelse(profs3$intensity > 0, profs3$intensity, profs3$mNZ)
profs3$logIntensity <- log10(profs3$intensity)

# Super LM:
AE_TC <- data.frame(profID = unique(profs3$profID), lrtOverall = NA, slope1 = NA, slope1p = NA,
                    slope34 = NA, slope34p = NA, AE34minAE1 = NA, AE34minAE1p = NA)
for(i in 1:nrow(AE_TC)){
  temp1 <- profs3 %>% filter(profID == AE_TC$profID[i])
  lm0 <- lm(logIntensity ~ 1, data = temp1)
  lm1 <- lm(logIntensity ~ whichSP * cycle, data = temp1)
  AE_TC$lrtOverall[i] <- anova(lm1, lm0)$`Pr(>F)`[2]
  # Coefficients for AE-1
  temp2 <- summary(multcomp::glht(lm1, linfct = matrix(c(0, 0, 1, 0), nrow = 1, byrow = TRUE)), 
          test = multcomp::adjusted("none"))$test
  AE_TC$slope1[i] <- temp2$coef
  AE_TC$slope1p[i] <- temp2$pvalues
  # Coefficients for AE-3/4:
  temp3 <- summary(multcomp::glht(lm1, linfct = matrix(c(0, 0, 1, 1), nrow = 1, byrow = TRUE)), 
                   test = multcomp::adjusted("none"))$test
  AE_TC$slope34[i] <- temp3$coef
  AE_TC$slope34p[i] <- temp3$pvalues
  AE_TC$AE34minAE1[i] <- - as.numeric(summary(emmeans::emmeans(lm1, pairwise ~ whichSP)$contrasts)["estimate"])
  AE_TC$AE34minAE1p[i] <- as.numeric(summary(emmeans::emmeans(lm1, pairwise ~ whichSP)$contrasts)["p.value"])
  print(i)
}

# Add profile m/z and RT:
AE_TC <- AE_TC %>% left_join(profInfo %>% 
                               select(profile_ID, inBlind = `in_blind?`, mean_mz, mean_RT), 
                             by = c("profID" = "profile_ID"))

# Adjusted p-values:
AE_TC$AE34minAE1q <- p.adjust(AE_TC$AE34minAE1p, method = "fdr")
AE_TC$slope1q <- p.adjust(AE_TC$slope1p, method = "fdr")
AE_TC$slope34q <- p.adjust(AE_TC$slope34p, method = "fdr")

# Join possible annotation data:
AE_TC$posFenIn <- ""
for(i in 1:nrow(AE_TC)){
  pMatch1 <- fenIn2$precursorMZ > AE_TC$mean_mz[i] - .001 & fenIn2$precursorMZ < AE_TC$mean_mz[i] + .001
  pMatch2 <- fenIn2[pMatch1, ]
  if(nrow(pMatch2) > 0){
    pMatch3 <- pMatch2$RT > AE_TC$mean_RT[i] / 60 - .5 & pMatch2$RT < AE_TC$mean_RT[i] / 60 + .5
    if(nrow(pMatch2[pMatch3,]) > 0){
      AE_TC$posFenIn[i] <- paste(pMatch2$feature[pMatch3], collapse = ";")
    }
  }
}

# Volcano plots:
# png(filename = paste0("./Plots/AE1_TC_Volcano_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 7, width = 8, units = "in", res = 600)
EnhancedVolcano::EnhancedVolcano(AE_TC, 
           lab = AE_TC$profID, x = "slope1", y = "slope1q", pCutoff = .05, 
           FCcutoff = 0.05, xlab = "Slope", ylab = expression(paste(-log[10],"(Adjusted p-value)")),
           legendPosition = "none", caption = "", pointSize = 1.25, labSize = 1.25, title = "", subtitle = "",
           xlim = c(-.15, .15), ylim = c(0, 3.5)) 
# dev.off()
# png(filename = paste0("./Plots/AE1_TC_VolcanoNOLab_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 7, width = 8, units = "in", res = 600)
EnhancedVolcano::EnhancedVolcano(AE_TC, 
     lab = "", x = "slope1", y = "slope1q", pCutoff = .05, 
     FCcutoff = 0.05, xlab = "Slope", ylab = expression(paste(-log[10],"(Adjusted p-value)")),
     legendPosition = "none", caption = "", title = "", subtitle = "",
     xlim = c(-.15, .15), ylim = c(0, 3.5)) 
# dev.off()

# png(filename = paste0("./Plots/AE34minusAE1_Volcano_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 7, width = 8, units = "in", res = 600)
EnhancedVolcano::EnhancedVolcano(AE_TC, 
     lab = AE_TC$profID, x = "AE34minAE1", y = "AE34minAE1q", pCutoff = .05, 
     FCcutoff = 0.5, xlab = "AE3/4 - AE1", ylab = expression(paste(-log[10],"(Adjusted p-value)")), 
     legendPosition = "none", caption = "",  ylim = c(0, 6.5), xlim = c(-2.5, 2.5), 
     pointSize = .95, labSize = .95, title = "", subtitle = "")
# dev.off()
# png(filename = paste0("./Plots/AE34minusAE1_VolcanoNOLab_",gsub("-", "", Sys.Date()), ".png"), 
#     height = 7, width = 8, units = "in", res = 600)
EnhancedVolcano::EnhancedVolcano(AE_TC, 
     lab = "", x = "AE34minAE1", y = "AE34minAE1q", pCutoff = .05, 
     FCcutoff = 0.5, xlab = "AE3/4 - AE1", ylab = expression(paste(-log[10],"(Adjusted p-value)")), 
     legendPosition = "none", caption = "", ylim = c(0, 6.5), xlim = c(-2.5, 2.5), 
     title = "", subtitle = "")
# dev.off()

# Specific profiles:
# png(filename = "./Plots/AE1_TC_13176.png", height = 4, width = 5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 13176 & whichSP == "AE-1"), aes(x = cycle, y = intensity)) + geom_point() + 
  stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #13176: AE-1 Only", subtitle = "m/z: 256.19064 (+/-0.17 ppm); RT: 5.1 min",
       x = "Cycle", y = "Intensity")
# dev.off()
# png(filename = "./Plots/AE1_TC_13176_Both.png", height = 4, width = 5.5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 13176), aes(x = cycle, color = whichSP, y = intensity)) + 
  geom_point() + stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #13176: AE-1 and AE-3/4", subtitle = "m/z: 256.19064 (+/-0.17 ppm); RT: 5.1 min",
       x = "Cycle", y = "Intensity", color = "Sampling\nPoint")
# dev.off()
# png(filename = "./Plots/AE1_TC_13176_Both2.png", height = 4, width = 5.5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 13176), aes(x = cycle, color = whichSP, y = log10(intensity))) + 
  geom_point() + stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #13176: AE-1 and AE-3/4 (Log scale)", subtitle = "m/z: 256.19064 (+/-0.17 ppm); RT: 5.1 min",
       x = "Cycle", y = expression(paste(log[10],"(Intensity)")), color = "Sampling\nPoint")
# dev.off()

ggplot(profs3 %>% filter(profID == 8845), aes(x = cycle, y = logIntensity, color = whichSP)) + geom_point() + 
  stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #8845", subtitle = "m/z: ; RT:  min",
       x = "Cycle", y = expression(paste(log[10],"(Intensity)")))

# Export time-course data:
writexl::write_xlsx(AE_TC, path = paste0("Results/AE_Timecourse_", gsub("-", "", Sys.Date()), ".xlsx"))

rm(lm0, lm1, p1, p2, p3, p4, p5, pca1, pca1DF, pMatch2, temp1, temp2, temp3, colEntropies,
   iSTDCV, iSTDDF, iSTDs, medByFile, medNorm, mISTD, mISTD2, mISTD3, mISTD4, profs2b, 
   pMatch1, pMatch3, i, profs3)

# Save:
save.image("working_20200702b.RData")

########### Algal present versus absent ###########
load("working_20200702b.RData")

# First get data for the profiles for AE:
profs3 <- profs2 %>% filter(grepl("AE-", fileName))
profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- ifelse(grepl("AE-3|AE-4", profs3$fileName), "AE-3/4", "AE-1")

# Add indicator for absent, intermediate, present:
profs3 <- profs3 %>% mutate(pres = case_when(intensity == 0 ~ "Absent", intensity > 5*10^6 ~ "Present", 
                                             TRUE ~ "Intermediate"))
profs3$pres <- factor(profs3$pres, levels = c("Absent", "Intermediate", "Present"), ordered = TRUE)

# Remove those absent in all:
# Starting count:
profs3 %>% select(profID) %>% unique() %>% nrow() # 48,255
profs3 <- profs3 %>% group_by(profID) %>% mutate(countPresent = sum(pres == "Present"), 
           removeThis = countPresent == 0) %>% filter(!removeThis) %>% select(-countPresent, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 12,236

AE_Pres <- data.frame(profID = unique(profs3$profID), lrtOverall = NA)
for(i in 1:nrow(AE_Pres)){
  temp1 <- profs3 %>% filter(profID == AE_Pres$profID[i])
  ordLogistic0 <- MASS::polr(pres ~ 1, data = temp1)
  ordLogistic1 <- MASS::polr(pres ~ whichSP, data = temp1)
  AE_Pres$lrtOverall[i] <- anova(ordLogistic0, ordLogistic1)$`Pr(Chi)`[2]
  print(i)
}

# Adjusted p-values:
AE_Pres$lrtOverallq <- p.adjust(AE_Pres$lrtOverall, method = "fdr")

# Add profile m/z and RT:
AE_Pres <- AE_Pres %>% left_join(profInfo %>% 
                               select(profile_ID, inBlind = `in_blind?`, mean_mz, mean_RT), 
                             by = c("profID" = "profile_ID"))

# Join possible annotation data:
AE_Pres$posFenIn <- ""
for(i in 1:nrow(AE_Pres)){
  pMatch1 <- fenIn2$precursorMZ > AE_Pres$mean_mz[i] - .001 & fenIn2$precursorMZ < AE_Pres$mean_mz[i] + .001
  pMatch2 <- fenIn2[pMatch1, ]
  if(nrow(pMatch2) > 0){
    pMatch3 <- pMatch2$RT > AE_Pres$mean_RT[i] / 60 - .5 & pMatch2$RT < AE_Pres$mean_RT[i] / 60 + .5
    if(nrow(pMatch2[pMatch3,]) > 0){
      AE_Pres$posFenIn[i] <- paste(pMatch2$feature[pMatch3], collapse = ";")
    }
  }
}

png(filename = "./Plots/AE_Present_4988.png", height = 5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(profs3 %>% filter(profID == 4988), aes(x = whichSP, color = whichSP, y = intensity, label = fileName)) + 
  geom_point() + geom_text_repel() + theme_bw() + 
  labs(title = "Profile #4988: AE-1 and AE-3/4", subtitle = "m/z: 182.0093; RT: 9.34 min",
       x = "Sampling Point", y = "Intensity", color = "Sampling\nPoint")
dev.off()
