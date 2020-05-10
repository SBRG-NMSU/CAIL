############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
oldPar <- par()

library(tidyverse)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/WaterProject"))

# Read in sample annotation:
sampleAnno <- read.csv("SampleAnno_20200509.csv")

# Load profile data:
load("enviMassOutput_20200509.RData")

# Set rownames on that data:
rownames(profs) <- sampleAnno$Name[match(rownames(profs), sampleAnno$ID)]

# Which are 0:
which(profs == 0, arr.ind = TRUE)

# Imputation for missing values:
impFun <- function(x){
  x[x == 0] <- min(x[x > 0]) / 2
  return(x)
}
profs <- apply(profs, 2, impFun)

# Convert to log-scale:
profs <- log2(profs)
profs <- scale(profs, center = TRUE, scale = FALSE)

# Entropy filter:
entropyFun <- function(x) entropy::entropy(entropy::discretize(x, 5), unit = "log2")
colEntropies <- apply(profs, 2, entropyFun)
hist(colEntropies)
profs <- profs[, colEntropies > 1]

heatmap(profs)

pca1 <- prcomp(profs)
plot(pca1)
pca1DF <- as.data.frame(pca1$x[, 1:2])
pca1DF$sampleID <- rownames(pca1DF)
pca1DF$sampleID <- gsub("90percent_", "", gsub("90per_", "", pca1DF$sampleID))
set.seed(3)
ggplot(pca1DF, aes(x = PC1, y = PC2, label = sampleID)) + geom_point() + geom_text_repel() + 
  theme_bw() + labs(title = "PCA prior to normalization")

########### IS-based normalization ###########
