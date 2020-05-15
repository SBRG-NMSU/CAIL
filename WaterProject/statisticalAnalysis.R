############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
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

# Imputation Function
impFun <- function(x){
  if(sum(x == 0) > 0){
    x[x == 0] <- min(x[x > 0]) / 2
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

############ Viz prior to normalization ############
# Heatmap
heatmap(aSProfs)

# PCA
pca1 <- prcomp(aSProfs)
plot(pca1)
pca1DF <- as.data.frame(pca1$x[, 1:2])
pca1DF$sampleID <- rownames(pca1DF)
pca1DF$sampleID <- gsub("90percent_", "", gsub("90per_", "", pca1DF$sampleID))
set.seed(3)
ggplot(pca1DF, aes(x = PC1, y = PC2, label = sampleID)) + geom_point() + geom_text_repel() + 
  theme_bw() + labs(title = "PCA prior to normalization")

########### Intensity distributions ###########
profs2 <- profs
profs2 <- profs2[, !apply(profs2, 2, misCheckFun2)]
profs2 <- apply(profs2, 2, impFun)
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName, levels = levels(sampleAnno$Name))
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
p1 <- ggplot(profs2, aes(x = fileName, y = log10(intensity))) + geom_boxplot() + 
  geom_hline(yintercept = median(log10(profs2$intensity)), color = "darkblue", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(x = "")

profs3 <- as.data.frame(profs)
profs3$fileName <- rownames(profs3)
profs3$fileName <- factor(profs3$fileName, levels = levels(sampleAnno$Name))
profs3 <- profs3 %>% gather(key = "profID", value = "intensity", -fileName)
profs3 <- profs3 %>% mutate(isFound = ifelse(intensity > 0, 1, 0))
profs3 <- profs3 %>% group_by(fileName) %>% summarize(profiles = sum(isFound))

p2 <- ggplot(profs3, aes(x = fileName, y = profiles)) + geom_point() + 
  geom_hline(yintercept = median(profs3$profiles), color = "darkred", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(x = "")

gridExtra::grid.arrange(p1, p2, nrow = 2)

########### IS-based normalization ###########
mISTD2 <- mISTD %>% group_by(file_ID, Name) %>% select(file_ID, Name, Intensity) %>% as.data.frame()
mISTD2 <- mISTD2 %>% left_join(sampleAnno %>% select(ID, fileName = Name, grp = tag3), by = c("file_ID" = "ID"))
iSTDCV <- mISTD2 %>% group_by(Name, grp) %>% summarize(cv = sd(Intensity) / mean(Intensity))

# Filter for the four that work (and remove fake QC samples):
mISTD3 <- mISTD2 %>% filter(Name %in% c(c("Acetaminophen-D4","Atenolol-D7", "Atrazine-D5", "Caffeine-13C3")))
mISTD3 <- mISTD3 %>% filter(grp != "FALSE")
medISTDInt <- mISTD3 %>% group_by(Name) %>% summarize(Mean = mean(Intensity, na.rm = TRUE))
mISTD3 <- mISTD3 %>% left_join(medISTDInt)
mISTD3$relDiff <- mISTD3$Intensity / mISTD3$Mean
mISTD4 <- mISTD3 %>% group_by(fileName) %>% summarize(invScale = mean(log10(relDiff)))

p4 <- ggplot(mISTD3, aes(x = Name, y = log10(Intensity), color = grp, group = fileName)) + geom_point() + 
  geom_line() + theme_bw()

# Plot of scaling:
ggplot(mISTD4, aes(x = fileName, y = invScale)) + geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + labs(x = "") 
mISTD4$invScale[mISTD4$fileName == "RO_3_2hr"] <- -.4
ggplot(mISTD4, aes(x = fileName, y = invScale)) + geom_point() + theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) + labs(x = "")

# Scale the data:
mISTD4$multFactor <- 10^(-mISTD4$invScale*.5)
sProfs <- profs
for(i in 1:nrow(sProfs)){
  sProfs[i, ] <- sProfs[i, ] * mISTD4$multFactor[match(rownames(sProfs)[i], mISTD4$fileName)]
}

# Scale the internal standards:
mISTD3 <- mISTD3 %>% left_join(mISTD4, by = "fileName") %>% mutate(Intensity2 = Intensity * multFactor)
p5 <- ggplot(mISTD3, aes(x = Name, y = log10(Intensity2), color = grp, group = fileName)) + geom_point() + 
  geom_line() + theme_bw()
gridExtra::grid.arrange(p4, p5, nrow = 2)

# Make a new boxplot
profs2 <- sProfs
profs2 <- profs2[, !apply(profs2, 2, misCheckFun2)]
profs2 <- apply(profs2, 2, impFun)
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName, levels = levels(sampleAnno$Name))
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
p3 <- ggplot(profs2, aes(x = fileName, y = log10(intensity))) + geom_boxplot() + 
  geom_hline(yintercept = median(log10(profs2$intensity)), color = "darkblue", lwd = 1) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90)) + labs(x = "")

gridExtra::grid.arrange(p1, p3, nrow = 2)

########### Algal time course ###########
profs2 <- sProfs[grepl("AE-1", rownames(sProfs)),]
profs2 <- profs2[, !apply(profs2, 2, misCheckFun2)]
profs2 <- apply(profs2, 2, impFun)
profs2 <- as.data.frame(profs2)
profs2$fileName <- rownames(profs2)
profs2$fileName <- factor(profs2$fileName, levels = levels(sampleAnno$Name))
profs2 <- profs2 %>% gather(key = "profID", value = "intensity", -fileName)
profs2$cycle <- as.integer(str_split(profs2$fileName, "_", simplify = TRUE)[,2])
profs2$logIntensity <- log10(profs2$intensity)

AE1_TC <- data.frame(profID = unique(profs2$profID), pValue = NA)
for(i in 1:nrow(AE1_TC)){
  lm1 <- lm(logIntensity ~ cycle, data = profs2 %>% filter(profID == AE1_TC$profID[i]))
  coefsLm1 <- coef(summary(lm1))
  AE1_TC$pValue[i] <- coefsLm1["cycle","Pr(>|t|)"]
}
ggplot(profs2 %>% filter(profID == 13176), aes(x = cycle, y = intensity)) + geom_point() + 
  geom_line() + theme_bw()
