############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
oldPar <- par()

library(tidyverse)
library(ggrepel)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/WaterProject"))

# Read in sample annotation:
sampleAnno <- read.csv("SampleAnno_20200509.csv")

# Load profile data:
load("enviMassOutput_20200509.RData")

# Set rownames on that data:
rownames(profs) <- sampleAnno$Name[match(rownames(profs), sampleAnno$ID)]

# Create data.frame with no blanks or "H20+IS samples":
actSamps <- sampleAnno$Name[sampleAnno$Type == "sample" & sampleAnno$profiled == TRUE]
aSProfs <- profs[rownames(profs) %in% actSamps,]

# Which are 0:
which(aSProfs == 0, arr.ind = TRUE)

# Remove all missing:
misCheckFun <- function(x){
  if(sum(x == 0) == length(x)){
    logical1 <- TRUE
  }else{
    logical1 <- FALSE
  }
  return(logical1)
}
aSProfs <- aSProfs[, !apply(aSProfs, 2, misCheckFun)]

# Imputation for missing values:
impFun <- function(x){
  if(sum(x == 0) > 0){
    x[x == 0] <- min(x[x > 0]) / 2
  }
  return(x)
}
aSProfs <- apply(aSProfs, 2, impFun)

# Convert to log-scale:
aSProfs <- log2(aSProfs)
aSProfs <- scale(aSProfs, center = TRUE, scale = FALSE)

# Entropy filter:
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, 5), unit = "log2")
colEntropies <- apply(aSProfs, 2, entropyFun)
hist(colEntropies)
aSProfs <- aSProfs[, colEntropies > 1]

heatmap(aSProfs)

pca1 <- prcomp(aSProfs)
plot(pca1)
pca1DF <- as.data.frame(pca1$x[, 1:2])
pca1DF$sampleID <- rownames(pca1DF)
pca1DF$sampleID <- gsub("90percent_", "", gsub("90per_", "", pca1DF$sampleID))
set.seed(3)
ggplot(pca1DF, aes(x = PC1, y = PC2, label = sampleID)) + geom_point() + geom_text_repel() + 
  theme_bw() + labs(title = "PCA prior to normalization")

########### IS-based normalization ###########
