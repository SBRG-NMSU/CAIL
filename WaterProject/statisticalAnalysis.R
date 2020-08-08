############ Prereqs ############
# Begin always run
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)
library(ggrepel)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/WaterProject"))
# End always run

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
sampleAnno <- read.csv("RData/SampleAnno_20200714.csv")
sampleAnno$ID <- as.character(sampleAnno$ID)
sampleAnno$Name <- factor(sampleAnno$Name, levels = sampleAnno$Name[sampleAnno$myOrder])

# Load profile data:
load("RData/enviMassOutput_20200714.RData")
rm(profs)
profs <- read.csv("RData/profiled_peaks.csv")
profInfo2 <- profs[, !grepl("X", names(profs))]
profs <- profs[, grepl("X", names(profs))]
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

# Import ID information:
fenIn <- readxl::read_excel("AE_MS2_CmpID/AE_CmpdID_20200620.xlsx", sheet = "FormulaEn_inSilico")
fenIn$feature <- as.character(as.integer(fenIn$feature))
# Unique:
fenIn2 <- fenIn %>% select(feature, precursorMZ, RT) %>% unique()

# Import MetFrag results:
load(file = "RData/metFragCand_20200804.RData")

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

# Remove pooled samples:
profs <- profs[!grepl("Pool", rownames(profs)),]

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

# png(filename = paste0("./Plots/PCA_prior2Norm_", gsub("-", "", Sys.Date()), ".png"),
#     height = 5.5, width = 8, units = "in", res = 600)
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
mISTD3 <- mISTD3 %>% filter(grp != "")
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

# png(filename = paste0("./Plots/iSTDIntens2_",gsub("-", "", Sys.Date()), ".png"),
#     height = 9, width = 7, units = "in", res = 600)
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

save.image("RData/working_20200714.RData")

########### Overall Present vs Absent ###########
load("RData/working_20200714.RData")

pvsA <- profs2 %>% left_join(sampleAnno %>% select(Name, ST = tag3), by = c("fileName" = "Name"))
pvsA$ST <- str_split(pvsA$ST, "-|_", simplify = TRUE)[,1]
pvsA <- pvsA %>% mutate(pres = case_when(intensity == 0 ~ "Absent", intensity > 5*10^6 ~ "Present", 
                                   TRUE ~ "Intermediate"))
pvsA <- pvsA %>% group_by(profID, ST, pres) %>% summarize(n = n())
pvsA <- pvsA %>% spread(key = pres, value = n)
pvsA$Absent[is.na(pvsA$Absent)] <- 0
pvsA$Intermediate[is.na(pvsA$Intermediate)] <- 0
pvsA$Present[is.na(pvsA$Present)] <- 0
pvsA$n <- pvsA$Absent + pvsA$Intermediate + pvsA$Present
pvsA$AllPresent <- pvsA$Present > .6 * pvsA$n

# 
pvsA2 <- pvsA %>% select(profID, ST, AllPresent) %>% group_by(ST, AllPresent) %>% summarize(n = sum(AllPresent)) %>%
  filter(AllPresent == TRUE) %>% select(-AllPresent)
pvsA %>% select(profID, ST, AllPresent) %>% group_by(ST, AllPresent) %>%
  filter(AllPresent == TRUE) %>% ungroup() %>% select(-AllPresent)

########### Algal time course ###########
load("RData/working_20200714.RData")

# First get data for the profiles for AE:
profs3 <- profs2 %>% filter(grepl("AE-", fileName))
profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- ifelse(grepl("AE-3|AE-4", profs3$fileName), "AE-3/4", "AE-1")

# Processing for time-course:
# 1. Remove profiles not observed in >40% of AE-1, AE-3, and AE-4 samples
# 2. Remove profiles with average peak area less than 5×10^6 in both AE-1 and AE-3/4

# Filter out super missing:
profs3 %>% select(profID) %>% unique() %>% nrow() # 53,672 profiles
profs3 <- profs3 %>% group_by(profID) %>% mutate(countNonZero = sum(intensity > 0), countTotal = n(),
                                                 removeIt = ifelse(countNonZero / countTotal <= .6, TRUE, FALSE)) %>%
  ungroup() %>% filter(!removeIt) %>% select(-(countNonZero:removeIt))
profs3 %>% select(profID) %>% unique() %>% nrow() # 1,206 profiles

# Filter out too low abundance: 
profs3 <- profs3 %>% group_by(profID, whichSP) %>% mutate(meanIntByGroup = mean(intensity), 
     lowAbundanceByGroup = meanIntByGroup < 5*10^6, groupN = n()) %>% ungroup() %>% group_by(profID) %>% 
     mutate(countLowAbundanceByGroup = sum(lowAbundanceByGroup), tooLow = countLowAbundanceByGroup == groupN) %>%
  group_by(profID) %>% mutate(removeThis = sum(tooLow) > 0) %>% ungroup() %>% filter(!removeThis) %>%
  select(-meanIntByGroup, -lowAbundanceByGroup, -groupN, -countLowAbundanceByGroup, -tooLow, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 760 profiles

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

# Adjusted p-values:
AE_TC$AE34minAE1q <- p.adjust(AE_TC$AE34minAE1p, method = "fdr")
AE_TC$slope1q <- p.adjust(AE_TC$slope1p, method = "fdr")
AE_TC$slope34q <- p.adjust(AE_TC$slope34p, method = "fdr")

# Add profile m/z and RT:
AE_TC <- AE_TC %>% left_join(profInfo2 %>% 
      select(profile_ID, profile_mean_mass, profile_mean_RT_min, homologue, neutral_mass, adduct), 
      by = c("profID" = "profile_ID")) %>% left_join(topHit, by = c("profID" = "prof")) %>% left_join(linksDF)

# Volcano plots:
# png(filename = paste0("./Plots/AE1_TC_Volcano_",gsub("-", "", Sys.Date()), ".png"),
#     height = 7, width = 8, units = "in", res = 600)
AE_TC2 <- AE_TC
AE_TC2$lab <- AE_TC2$profID
AE_TC2$lab[-log10(AE_TC2$slope1q) < 1.30103 | abs(AE_TC2$slope1) < .05] <- ""
set.seed(3)
ggplot(AE_TC2, aes(x = slope1, y = -log10(slope1q), label = lab)) + 
  geom_point(pch = 21,color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept = 1.30103, lty = 2) + geom_vline(xintercept = -.05, lty = 2) + 
  geom_vline(xintercept = .05, lty = 2) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5, segment.size = .5) +
  theme_bw() + labs(x = "AE-1 Time-Course Slope", y = "-Log10(q-value)")
# dev.off()

# png(filename = paste0("./Plots/AE34minusAE1_Volcano_",gsub("-", "", Sys.Date()), ".png"),
#     height = 7, width = 8, units = "in", res = 600)
AE_TC2 <- AE_TC
AE_TC2$lab <- AE_TC2$profID
AE_TC2$lab[-log10(AE_TC2$AE34minAE1q) < 1.30103 | abs(AE_TC2$AE34minAE1) < .75] <- ""
set.seed(3)
ggplot(AE_TC2, aes(x = AE34minAE1, y = -log10(AE34minAE1q), label = lab)) + 
  geom_point(pch = 21,color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept = 1.30103, lty = 2) + geom_vline(xintercept = -.75, lty = 2) + 
  geom_vline(xintercept = .75, lty = 2) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5, segment.size = .5) +
  theme_bw() + xlim(-2, 2.25) + labs(x = "AE-3/4 - AE-1", y = "-Log10(q-value)")
# dev.off()

# Specific profiles:
# png(filename = "./Plots/AE1_TC_16791.png", height = 4, width = 5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 16791 & whichSP == "AE-1"), aes(x = cycle, y = intensity)) + geom_point() + 
  stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #16791: AE-1 Only", subtitle = "m/z: 256.19064 (+/-0.17 ppm); RT: 5.1 min",
       x = "Cycle", y = "Intensity")
# dev.off()
# png(filename = "./Plots/AE1_TC_16791_Both.png", height = 4, width = 5.5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 16791), aes(x = cycle, color = whichSP, y = intensity)) + 
  geom_point() + stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #16791: AE-1 and AE-3/4", subtitle = "m/z: 256.19064 (+/-0.17 ppm); RT: 5.1 min",
       x = "Cycle", y = "Intensity", color = "Sampling\nPoint")
# dev.off()
# png(filename = "./Plots/AE1_TC_16791_Both2.png", height = 4, width = 5.5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 16791), aes(x = cycle, color = whichSP, y = log10(intensity))) + 
  geom_point() + stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #16791: AE-1 and AE-3/4 (Log scale)", subtitle = "m/z: 256.19064 (+/-0.17 ppm); RT: 5.1 min",
       x = "Cycle", y = expression(paste(log[10],"(Intensity)")), color = "Sampling\nPoint")
# dev.off()

# Export time-course data:
AE_TC <- AE_TC %>% select(profID, mz = profile_mean_mass, rt = profile_mean_RT_min, neutral_mass,
                          adduct, formula, howID = id, topCandidate = name, links,
                          slope1, slope34, AE34minusAE1 = AE34minAE1, slope1p, slope34p, AE34minusAE1p = AE34minAE1p,
                          slope1q, slope34q, AE34minusAE1q = AE34minAE1q)
writexl::write_xlsx(AE_TC, path = paste0("Results/AE_Timecourse_", gsub("-", "", Sys.Date()), ".xlsx"))

rm(lm0, lm1, p1, p2, p3, p4, p5, pca1, pca1DF, pMatch2, temp1, temp2, temp3, colEntropies,
   iSTDCV, iSTDDF, iSTDs, medByFile, medNorm, mISTD, mISTD2, mISTD3, mISTD4, profs2b, 
   pMatch1, pMatch3, i, profs3, AE_TC2)

# Save:
save.image("RData/working_20200714b.RData")

########### Algal present versus absent ###########
load("RData/working_20200714b.RData")

# First get data for the profiles for AE:
profs3 <- profs2 %>% filter(grepl("AE-", fileName))
profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- factor(ifelse(grepl("AE-3|AE-4", profs3$fileName), "AE-3/4", "AE-1"))

# Add indicator for absent, intermediate, present:
profs3 <- profs3 %>% mutate(pres = case_when(intensity == 0 ~ "Absent", intensity > 5*10^6 ~ "Present", 
                                             TRUE ~ "Intermediate"))
profs3$pres <- factor(profs3$pres, levels = c("Absent", "Intermediate", "Present"), ordered = TRUE)

# Remove those absent in all:
# Starting count:
profs3 %>% select(profID) %>% unique() %>% nrow() # 53,672
profs3 <- profs3 %>% group_by(profID) %>% mutate(countPresent = sum(pres == "Present"), 
           removeThis = countPresent == 0) %>% filter(!removeThis) %>% select(-countPresent, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 7,850

Pres <- data.frame(profID = unique(profs3$profID), diff = NA, pValue = NA)
Pres2 <- list()
for(i in 1:nrow(Pres)){
  # Find all intensity for an individual feature:
  temp1 <- profs3 %>% filter(profID == Pres$profID[i])
  
  # Tabulate present vs. absent:
  temp2 <- as.data.frame(xtabs(~pres + whichSP, data = temp1)) 
  temp2$pres <- paste0(temp2$pres, "_", temp2$whichSP)
  temp2 <- temp2 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp2) <- paste0(names(temp2), "_Freq")
  
  # Tabulate proportions:
  temp3 <- as.data.frame(prop.table(xtabs(~pres + whichSP, data = temp1), margin = 2))
  temp3$pres <- paste0(temp3$pres, "_", temp3$whichSP)
  temp3 <- temp3 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp3) <- paste0(names(temp3), "_Prop")
  
  # Independence test:
  if(length(table(as.integer(temp1$pres))) > 1){
    test1 <- coin::independence_test(pres ~ whichSP, data = temp1)
    Pres$pValue[i] <- pnorm(abs(test1@statistic@teststatistic), lower.tail = FALSE) * 2
    temp1$pres2 <- as.integer(temp1$pres)
    wc2 <- temp1 %>% group_by(whichSP) %>% summarize(mean = mean(pres2))
    Pres$diff[i] <- wc2$mean[2] - wc2$mean[1]
  }else{
    Pres$diff[i] <- 0
    Pres$pValue[i] <- 1
  }
  
  # Eports:
  Pres2[[i]] <- cbind(profID = Pres$profID[i], temp2, temp3)
  print(i)
}
Pres2 <- do.call("rbind", Pres2)
Pres <- Pres %>% left_join(Pres2)

# Adjusted p-values:
Pres$qValue <- p.adjust(Pres$pValue, method = "fdr")

# Join with annotation data:
Pres <- Pres %>% left_join(profInfo2 %>% 
   select(profile_ID, profile_mean_mass, profile_mean_RT_min, homologue, neutral_mass, adduct), 
   by = c("profID" = "profile_ID")) %>% left_join(topHit, by = c("profID" = "prof")) %>% left_join(linksDF)
Pres <- Pres %>% select(profID, mz = profile_mean_mass, rt = profile_mean_RT_min, neutral_mass,
                 adduct, formula, howID = id, topCandidate = name, links,
                `Absent_AE-1_Freq`, `Absent_AE-3/4_Freq`, `Intermediate_AE-1_Freq`, `Intermediate_AE-3/4_Freq`,
                `Present_AE-1_Freq`, `Present_AE-3/4_Freq`, `Absent_AE-1_Prop`, `Absent_AE-3/4_Prop`,
                `Intermediate_AE-1_Prop`, `Intermediate_AE-3/4_Prop`, `Present_AE-1_Prop`, `Present_AE-3/4_Prop`, 
                MeanDifference = diff, pValue, qValue)

# png(filename = paste0("./Plots/AEPresent_Volcano_",gsub("-", "", Sys.Date()), ".png"),
#     height = 7, width = 8, units = "in", res = 600)
Pres2 <- Pres
Pres2$lab <- Pres2$profID
Pres2$lab[-log10(Pres2$qValue) < -log10(.1) | abs(Pres2$MeanDifference) < 1.25] <- ""
set.seed(33)
ggplot(Pres2, aes(x = MeanDifference, y = -log10(qValue), label = lab)) + 
  geom_point(pch = 21,color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept =-log10(.1), lty = 2) + geom_vline(xintercept = -1.25, lty = 2) + 
  geom_vline(xintercept = 1.25, lty = 2) + 
  geom_text_repel(size = 1.5, segment.colour = "grey30", segment.alpha = .35, segment.size = .5) +
  theme_bw() + labs(x = "Mean Difference", y = "-Log10(q-value)")
# dev.off()

# png(filename = "./Plots/AE_Present_45881.png", height = 5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(profs3 %>% filter(profID == 45881), aes(x = whichSP, color = whichSP, y = intensity, label = fileName)) + 
  geom_point() + geom_text_repel() + theme_bw() + 
  labs(title = "Profile #45881: AE-1 and AE-3/4", subtitle = "m/z: 453.3436; RT: 7.09 min",
       x = "Sampling Point", y = "Intensity", color = "Sampling\nPoint")
# dev.off()

# Export:
writexl::write_xlsx(Pres, path = paste0("Results/AE_Pres_", gsub("-", "", Sys.Date()), ".xlsx"))

# Save and cleanup:
Pres_AE <- Pres
rm(Pres, Pres2, wc1, wc2, temp1, temp2, temp3, profs3, test1)
save.image("RData/working_20200702c.RData")

########### Secondary Effluent to RO ###########
load("RData/working_20200702c.RData")

# First get data for the profiles for Secondary Effluent and RO:
profs3 <- profs2 %>% filter((grepl("SE-", fileName) & !grepl("Pool", fileName)) | 
                              (grepl("RO_", fileName) & !grepl("Secondary", fileName)))
# Toss out bad RO samples:
profs3 <- profs3 %>% filter(!fileName %in% paste0("RO_", 9:14))
profs3$fileName <- factor(profs3$fileName)
table(profs3$fileName)

profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- ifelse(grepl("RO", profs3$fileName), "RO", "SE")
profs3$whichSP <- factor(profs3$whichSP)

# Add indicator for absent, intermediate, present:
profs3 <- profs3 %>% mutate(pres = case_when(intensity == 0 ~ "Absent", intensity > 5*10^6 ~ "Present", 
                                             TRUE ~ "Intermediate"))
profs3$pres <- factor(profs3$pres, levels = c("Absent", "Intermediate", "Present"), ordered = TRUE)

# Remove those absent in all:
# Starting count:
profs3 %>% select(profID) %>% unique() %>% nrow() # 53,672
profs3 <- profs3 %>% group_by(profID) %>% mutate(countPresent = sum(pres == "Present"), 
                                                 removeThis = countPresent == 0) %>% filter(!removeThis) %>% select(-countPresent, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 11,499

Pres <- data.frame(profID = unique(profs3$profID), diff = NA, pValue = NA)
Pres2 <- list()
for(i in 1:nrow(Pres)){
  # Find all intensity for an individual feature:
  temp1 <- profs3 %>% filter(profID == Pres$profID[i])
  
  # Tabulate present vs. absent:
  temp2 <- as.data.frame(xtabs(~pres + whichSP, data = temp1)) 
  temp2$pres <- paste0(temp2$pres, "_", temp2$whichSP)
  temp2 <- temp2 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp2) <- paste0(names(temp2), "_Freq")
  
  # Tabulate proportions:
  temp3 <- as.data.frame(prop.table(xtabs(~pres + whichSP, data = temp1), margin = 2))
  temp3$pres <- paste0(temp3$pres, "_", temp3$whichSP)
  temp3 <- temp3 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp3) <- paste0(names(temp3), "_Prop")
  
  # Independence test:
  if(length(table(as.integer(temp1$pres))) > 1){
    test1 <- coin::independence_test(pres ~ whichSP, data = temp1)
    Pres$pValue[i] <- pnorm(abs(test1@statistic@teststatistic), lower.tail = FALSE) * 2
    temp1$pres2 <- as.integer(temp1$pres)
    wc2 <- temp1 %>% group_by(whichSP) %>% summarize(mean = mean(pres2))
    Pres$diff[i] <- wc2$mean[2] - wc2$mean[1]
  }else{
    Pres$diff[i] <- 0
    Pres$pValue[i] <- 1
  }
  
  # Eports:
  Pres2[[i]] <- cbind(profID = Pres$profID[i], temp2, temp3)
  print(i)
}
Pres2 <- do.call("rbind", Pres2)
Pres <- Pres %>% left_join(Pres2)

# Adjusted p-values:
Pres$qValue <- p.adjust(Pres$pValue, method = "fdr")

# Join with annotation data:
Pres <- Pres %>% left_join(profInfo2 %>% 
     select(profile_ID, profile_mean_mass, profile_mean_RT_min, homologue, neutral_mass, adduct), 
     by = c("profID" = "profile_ID")) %>% left_join(topHit, by = c("profID" = "prof")) %>% left_join(linksDF)
Pres <- Pres %>% select(profID, mz = profile_mean_mass, rt = profile_mean_RT_min, neutral_mass,
        adduct, formula, howID = id, topCandidate = name, links,
        `Absent_SE_Freq`, `Absent_RO_Freq`, `Intermediate_SE_Freq`, `Intermediate_RO_Freq`,
        `Present_SE_Freq`, `Present_RO_Freq`, `Absent_SE_Prop`, `Absent_RO_Prop`,
        `Intermediate_SE_Prop`, `Intermediate_RO_Prop`, `Present_SE_Prop`, `Present_RO_Prop`, 
        MeanDifference = diff, pValue, qValue)

# png(filename = paste0("./Plots/SE_RO_Present_Volcano_",gsub("-", "", Sys.Date()), ".png"),
#     height = 7, width = 8, units = "in", res = 600)
Pres2 <- Pres
Pres2 <- Pres2 %>% group_by(MeanDifference, qValue) %>% mutate(lab = n())
Pres2$lab[-log10(Pres2$qValue) < -log10(0.001) | abs(Pres2$MeanDifference) < 1.5] <- ""
Pres2 <- Pres2 %>% select(MeanDifference, qValue, lab) %>% unique()
set.seed(33)
ggplot(Pres2, aes(x = MeanDifference, y = -log10(qValue), label = lab)) + 
  geom_point(pch = 21, color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept =-log10(.001), lty = 2) + geom_vline(xintercept = -1.5, lty = 2) + 
  geom_vline(xintercept = 1.5, lty = 2) + 
  geom_text_repel(size = 3, segment.colour = "grey30", segment.alpha = .35, segment.size = .5) +
  theme_bw() + labs(x = "Mean Difference", y = "-Log10(q-value)")
# dev.off()

# png(filename = "./Plots/SE_RO_50460.png", height = 5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(profs3 %>% filter(profID == 50460), aes(x = whichSP, color = whichSP, y = intensity, label = fileName)) + 
  geom_point() + geom_text_repel(size = 1.25, segment.size = .25, segment.alpha = .25) + theme_bw() + 
  labs(title = "Profile #50460: SE and RO", subtitle = "m/z: 507.271942; RT:  10.56 min",
       x = "Sampling Point", y = "Intensity", color = "Sampling\nPoint")
# dev.off()

# Export:
writexl::write_xlsx(Pres, path = paste0("Results/SE_Product_Pres_", gsub("-", "", Sys.Date()), ".xlsx"))

# Save and cleanup:
Pres_SE_RO <- Pres
rm(Pres, Pres2, test1, temp1, temp2, temp3, profs3)
save.image("RData/working_20200702d.RData")

########### Secondary Effluent to Secondary RO ###########
load("RData/working_20200702d.RData")

# First get data for the profiles for Secondary Effluent and product water:
profs3 <- profs2 %>% filter((grepl("SE-", fileName) & !grepl("Pool", fileName)) | 
                              (grepl("RO_", fileName) & grepl("Secondary", fileName)))
profs3$fileName <- factor(profs3$fileName)
table(profs3$fileName)

profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- ifelse(grepl("SecondaryRO", profs3$fileName), "SecondaryRO", "SE")
profs3$whichSP <- factor(profs3$whichSP)

# Add indicator for absent, intermediate, present:
profs3 <- profs3 %>% mutate(pres = case_when(intensity == 0 ~ "Absent", intensity > 5*10^6 ~ "Present", 
                                             TRUE ~ "Intermediate"))
profs3$pres <- factor(profs3$pres, levels = c("Absent", "Intermediate", "Present"), ordered = TRUE)

# Remove those absent in all:
# Starting count:
profs3 %>% select(profID) %>% unique() %>% nrow() # 53,672
profs3 <- profs3 %>% group_by(profID) %>% mutate(countPresent = sum(pres == "Present"), 
                                                 removeThis = countPresent == 0) %>% filter(!removeThis) %>% select(-countPresent, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 10,574

Pres <- data.frame(profID = unique(profs3$profID), diff = NA, pValue = NA)
Pres2 <- list()
for(i in 1:nrow(Pres)){
  # Find all intensity for an individual feature:
  temp1 <- profs3 %>% filter(profID == Pres$profID[i])
  
  # Tabulate present vs. absent:
  temp2 <- as.data.frame(xtabs(~pres + whichSP, data = temp1)) 
  temp2$pres <- paste0(temp2$pres, "_", temp2$whichSP)
  temp2 <- temp2 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp2) <- paste0(names(temp2), "_Freq")
  
  # Tabulate proportions:
  temp3 <- as.data.frame(prop.table(xtabs(~pres + whichSP, data = temp1), margin = 2))
  temp3$pres <- paste0(temp3$pres, "_", temp3$whichSP)
  temp3 <- temp3 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp3) <- paste0(names(temp3), "_Prop")
  
  # Independence test:
  if(length(table(as.integer(temp1$pres))) > 1){
    test1 <- coin::independence_test(pres ~ whichSP, data = temp1)
    Pres$pValue[i] <- pnorm(abs(test1@statistic@teststatistic), lower.tail = FALSE) * 2
    temp1$pres2 <- as.integer(temp1$pres)
    wc2 <- temp1 %>% group_by(whichSP) %>% summarize(mean = mean(pres2))
    Pres$diff[i] <- wc2$mean[2] - wc2$mean[1]
  }else{
    Pres$diff[i] <- 0
    Pres$pValue[i] <- 1
  }
  
  # Eports:
  Pres2[[i]] <- cbind(profID = Pres$profID[i], temp2, temp3)
  print(i)
}
Pres2 <- do.call("rbind", Pres2)
Pres <- Pres %>% left_join(Pres2)

# Adjusted p-values:
Pres$qValue <- p.adjust(Pres$pValue, method = "fdr")

# Join with annotation data:
Pres <- Pres %>% left_join(profInfo2 %>% 
                             select(profile_ID, profile_mean_mass, profile_mean_RT_min, homologue, neutral_mass, adduct), 
                           by = c("profID" = "profile_ID")) %>% left_join(topHit, by = c("profID" = "prof")) %>% left_join(linksDF)

Pres <- Pres %>% select(profID, mz = profile_mean_mass, rt = profile_mean_RT_min, neutral_mass,
                        adduct, formula, howID = id, topCandidate = name, links,
                        `Absent_SE_Freq`, `Absent_SecondaryRO_Freq`, `Intermediate_SE_Freq`, `Intermediate_SecondaryRO_Freq`,
                        `Present_SE_Freq`, `Present_SecondaryRO_Freq`, `Absent_SE_Prop`, `Absent_SecondaryRO_Prop`,
                        `Intermediate_SE_Prop`, `Intermediate_SecondaryRO_Prop`, `Present_SE_Prop`, `Present_SecondaryRO_Prop`, 
                        MeanDifference = diff, pValue, qValue)

# Volcano Plot:
# png(filename = paste0("./Plots/SE_SecondaryRO_Present_Volcano_",gsub("-", "", Sys.Date()), ".png"),
#     height = 7, width = 8, units = "in", res = 600)
Pres2 <- Pres
Pres2 <- Pres2 %>% group_by(MeanDifference, qValue) %>% mutate(lab = n())
Pres2$lab[-log10(Pres2$qValue) < -log10(0.001) | abs(Pres2$MeanDifference) < 1.5] <- ""
Pres2 <- Pres2 %>% select(MeanDifference, qValue, lab) %>% unique()
set.seed(33)
ggplot(Pres2, aes(x = MeanDifference, y = -log10(qValue), label = lab)) + 
  geom_point(pch = 21, color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept =-log10(.001), lty = 2) + geom_vline(xintercept = -1.5, lty = 2) + 
  geom_vline(xintercept = 1.5, lty = 2) + 
  geom_text_repel(size = 3, segment.colour = "grey30", segment.alpha = .35, segment.size = .5) +
  theme_bw() + labs(x = "Mean Difference", y = "-Log10(q-value)")
# dev.off()

# png(filename = "./Plots/SE_SecondaryRO_16427.png", height = 5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(profs3 %>% filter(profID == 16427), aes(x = whichSP, color = whichSP, y = intensity, label = fileName)) + 
  geom_point() + geom_text_repel(size = 1.25, segment.size = .25, segment.alpha = .25) + theme_bw() + 
  labs(title = "Profile #16427: SE and Secondary RO", subtitle = "m/z: 254.2478; RT: 10.69 min",
       x = "Sampling Point", y = "Intensity", color = "Sampling\nPoint")
# dev.off()

# Export:
writexl::write_xlsx(Pres, path = paste0("Results/SE_SecondaryRO_Pres_", gsub("-", "", Sys.Date()), ".xlsx"))

# Save and cleanup:
Pres_SE_SecondaryRO <- Pres
rm(Pres, Pres2, test1, temp1, temp2, temp3, profs3)
save.image("RData/working_20200702e.RData")

########### Algal Effluent to RO ###########
load("RData/working_20200702e.RData")

# First get data for the profiles for Secondary Effluent and product water:
profs3 <- profs2 %>% filter((grepl("AE-", fileName) & !grepl("Pool", fileName)) | 
                              (grepl("RO_", fileName) & !grepl("Secondary", fileName)))

# Toss out bad RO samples:
profs3 <- profs3 %>% filter(!fileName %in% paste0("RO_", 9:14))
profs3$fileName <- factor(profs3$fileName)
table(profs3$fileName)

profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)
profs3$whichSP <- ifelse(grepl("RO", profs3$fileName), "RO", "AE")
profs3$whichSP <- factor(profs3$whichSP)

# Add indicator for absent, intermediate, present:
profs3 <- profs3 %>% mutate(pres = case_when(intensity == 0 ~ "Absent", intensity > 5*10^6 ~ "Present", 
                                             TRUE ~ "Intermediate"))
profs3$pres <- factor(profs3$pres, levels = c("Absent", "Intermediate", "Present"), ordered = TRUE)

# Remove those absent in all:
# Starting count:
profs3 %>% select(profID) %>% unique() %>% nrow() # 53,672
profs3 <- profs3 %>% group_by(profID) %>% mutate(countPresent = sum(pres == "Present"), 
                                                 removeThis = countPresent == 0) %>% filter(!removeThis) %>% select(-countPresent, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 14,125

Pres <- data.frame(profID = unique(profs3$profID), diff = NA, pValue = NA)
Pres2 <- list()
for(i in 1:nrow(Pres)){
  # Find all intensity for an individual feature:
  temp1 <- profs3 %>% filter(profID == Pres$profID[i])
  
  # Tabulate present vs. absent:
  temp2 <- as.data.frame(xtabs(~pres + whichSP, data = temp1)) 
  temp2$pres <- paste0(temp2$pres, "_", temp2$whichSP)
  temp2 <- temp2 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp2) <- paste0(names(temp2), "_Freq")
  
  # Tabulate proportions:
  temp3 <- as.data.frame(prop.table(xtabs(~pres + whichSP, data = temp1), margin = 2))
  temp3$pres <- paste0(temp3$pres, "_", temp3$whichSP)
  temp3 <- temp3 %>% select(-whichSP) %>% spread(key = pres, value = Freq)
  names(temp3) <- paste0(names(temp3), "_Prop")
  
  # Independence test:
  if(length(table(as.integer(temp1$pres))) > 1){
    test1 <- coin::independence_test(pres ~ whichSP, data = temp1)
    Pres$pValue[i] <- pnorm(abs(test1@statistic@teststatistic), lower.tail = FALSE) * 2
    temp1$pres2 <- as.integer(temp1$pres)
    wc2 <- temp1 %>% group_by(whichSP) %>% summarize(mean = mean(pres2))
    Pres$diff[i] <- wc2$mean[2] - wc2$mean[1]
  }else{
    Pres$diff[i] <- 0
    Pres$pValue[i] <- 1
  }
  
  # Eports:
  Pres2[[i]] <- cbind(profID = Pres$profID[i], temp2, temp3)
  print(i)
}
Pres2 <- do.call("rbind", Pres2)
Pres <- Pres %>% left_join(Pres2)

# Adjusted p-values:
Pres$qValue <- p.adjust(Pres$pValue, method = "fdr")

# Join with annotation data:
Pres <- Pres %>% left_join(profInfo2 %>% 
                             select(profile_ID, profile_mean_mass, profile_mean_RT_min, homologue, neutral_mass, adduct), 
                           by = c("profID" = "profile_ID")) %>% left_join(topHit, by = c("profID" = "prof")) %>% left_join(linksDF)

Pres <- Pres %>% select(profID, mz = profile_mean_mass, rt = profile_mean_RT_min, neutral_mass,
                        adduct, formula, howID = id, topCandidate = name, links,
                        `Absent_AE_Freq`, `Absent_RO_Freq`, `Intermediate_AE_Freq`, `Intermediate_RO_Freq`,
                        `Present_AE_Freq`, `Present_RO_Freq`, `Absent_AE_Prop`, `Absent_RO_Prop`,
                        `Intermediate_AE_Prop`, `Intermediate_RO_Prop`, `Present_AE_Prop`, `Present_RO_Prop`, 
                        MeanDifference = diff, pValue, qValue)

# Volcano Plot:
# png(filename = paste0("./Plots/AE_RO_Present_Volcano_",gsub("-", "", Sys.Date()), ".png"),
#     height = 7, width = 8, units = "in", res = 600)
Pres2 <- Pres
Pres2 <- Pres2 %>% group_by(MeanDifference, qValue) %>% mutate(lab = n())
Pres2$lab[-log10(Pres2$qValue) < -log10(0.001) | abs(Pres2$MeanDifference) < 1.5] <- ""
Pres2 <- Pres2 %>% select(MeanDifference, qValue, lab) %>% unique()
set.seed(33)
ggplot(Pres2, aes(x = MeanDifference, y = -log10(qValue), label = lab)) + 
  geom_point(pch = 21, color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept =-log10(.001), lty = 2) + geom_vline(xintercept = -1.5, lty = 2) + 
  geom_vline(xintercept = 1.5, lty = 2) + 
  geom_text_repel(size = 3, segment.colour = "grey30", segment.alpha = .35, segment.size = .5) +
  theme_bw() + labs(x = "Mean Difference", y = "-Log10(q-value)")
# dev.off()

# png(filename = "./Plots/AE_RO_1992.png", height = 5, width = 6, units = "in", res = 600)
set.seed(3)
ggplot(profs3 %>% filter(profID == 1992), aes(x = whichSP, color = whichSP, y = intensity, label = fileName)) + 
  geom_point() + geom_text_repel(size = 1.25, segment.size = .25, segment.alpha = .25) + theme_bw() + 
  labs(title = "Profile #1992: AE and RO", subtitle = "m/z: 136.0732; RT: 5.66 min",
       x = "Sampling Point", y = "Intensity", color = "Sampling\nPoint")
# dev.off()

# Export:
writexl::write_xlsx(Pres, path = paste0("Results/AE_RO_Pres_", gsub("-", "", Sys.Date()), ".xlsx"))

# Save and cleanup:
Pres_AE_RO <- Pres
rm(Pres, Pres2, test1, temp1, temp2, temp3, profs3)
save.image("RData/working_20200702f.RData")

########### RO Time Course ###########
load("RData/working_20200702f.RData")

# First get data for the profiles for RO:
profs3 <- profs2 %>% filter(grepl("RO_", fileName) & !grepl("Secondary", fileName))
profs3 <- profs3 %>% filter(!fileName %in% paste0("RO_", 9:14))
profs3$fileName <- factor(profs3$fileName)
table(profs3$fileName)

profs3$cycle <- as.integer(str_split(profs3$fileName, "_", simplify = TRUE)[,2])
profs3$logIntensity <- log10(profs3$intensity)

# Processing for time-course:
# 1. Remove profiles not observed in >40% of samples
# 2. Remove profiles with average peak area less than 5×10^6 in samples

# Filter out super missing:
profs3 %>% select(profID) %>% unique() %>% nrow() # 53,672 profiles
profs3 <- profs3 %>% group_by(profID) %>% mutate(countNonZero = sum(intensity > 0), countTotal = n(),
                                                 removeIt = ifelse(countNonZero / countTotal <= .6, TRUE, FALSE)) %>%
  ungroup() %>% filter(!removeIt) %>% select(-(countNonZero:removeIt))
profs3 %>% select(profID) %>% unique() %>% nrow() # 2,783 profiles

# Filter out too low abundance: 
profs3 <- profs3 %>% group_by(profID) %>% mutate(meanIntByGroup = mean(intensity), 
                                                          lowAbundanceByGroup = meanIntByGroup < 5*10^6, groupN = n()) %>% ungroup() %>% group_by(profID) %>% 
  mutate(countLowAbundanceByGroup = sum(lowAbundanceByGroup), tooLow = countLowAbundanceByGroup == groupN) %>%
  group_by(profID) %>% mutate(removeThis = sum(tooLow) > 0) %>% ungroup() %>% filter(!removeThis) %>%
  select(-meanIntByGroup, -lowAbundanceByGroup, -groupN, -countLowAbundanceByGroup, -tooLow, -removeThis)
profs3 %>% select(profID) %>% unique() %>% nrow() # 2,091 profiles

# Imputation
profs3 <- profs3 %>% group_by(profID) %>% mutate(mNZ = minNonZero(intensity))
profs3$intensity <- ifelse(profs3$intensity > 0, profs3$intensity, profs3$mNZ)
profs3$logIntensity <- log10(profs3$intensity)

# Super LM:
TC <- data.frame(profID = unique(profs3$profID), lrtOverall = NA, slope = NA, slopeP = NA)
for(i in 1:nrow(TC)){
  temp1 <- profs3 %>% filter(profID == TC$profID[i])
  lm0 <- lm(logIntensity ~ 1, data = temp1)
  lm1 <- lm(logIntensity ~ cycle, data = temp1)
  TC$lrtOverall[i] <- anova(lm1, lm0)$`Pr(>F)`[2]
  # Coefficients for beta:
  temp2 <- as.data.frame(coef(summary(lm1)))
  
  TC$slope[i] <- temp2["cycle","Estimate"]
  TC$slopeP[i] <- temp2["cycle","Pr(>|t|)"]
  print(i)
}

# Adjusted p-values:
TC$slopeQ <- p.adjust(TC$slopeP, method = "fdr")

# Add profile m/z and RT:
TC <- TC %>% left_join(profInfo2 %>% 
        select(profile_ID, profile_mean_mass, profile_mean_RT_min, homologue, neutral_mass, adduct), 
        by = c("profID" = "profile_ID")) %>% left_join(topHit, by = c("profID" = "prof")) %>% left_join(linksDF)

# Volcano plots:
png(filename = paste0("./Plots/RO_TC_Volcano_",gsub("-", "", Sys.Date()), ".png"),
    height = 7, width = 8, units = "in", res = 600)
TC2 <- TC
TC2$lab <- TC2$profID
TC2$lab[-log10(TC2$slopeQ) < 2 | abs(TC2$slope) < .075] <- ""
set.seed(3)
ggplot(TC2, aes(x = slope, y = -log10(slopeQ), label = lab)) + 
  geom_point(pch = 21,color = "grey30", fill = "dodgerblue", alpha = .5) + 
  geom_hline(yintercept = 2, lty = 2) + geom_vline(xintercept = -.075, lty = 2) + 
  geom_vline(xintercept = .075, lty = 2) + 
  geom_text_repel(size = 2, segment.colour = "grey30", segment.alpha = .5, segment.size = .5) +
  theme_bw() + labs(x = "RO Time-Course Slope", y = "-Log10(q-value)")
dev.off()

# png(filename = "./Plots/RO_TC_42274.png", height = 4, width = 5, units = "in", res = 600)
ggplot(profs3 %>% filter(profID == 42274), aes(x = cycle, y = log10(intensity))) + geom_point() + 
  stat_smooth(method = "lm") + theme_bw() + 
  labs(title = "Profile #42274", subtitle = "m/z: 421.0366; RT: 7.47 min",
       x = "Cycle", y = "Log10(Intensity)")
# dev.off()

# Export time-course data:
TC <- TC %>% select(profID, mz = profile_mean_mass, rt = profile_mean_RT_min, neutral_mass,
                          adduct, formula, howID = id, topCandidate = name, links,
                          slope, slopeP, slopeQ)
writexl::write_xlsx(TC, path = paste0("Results/RO_Timecourse_", gsub("-", "", Sys.Date()), ".xlsx"))
