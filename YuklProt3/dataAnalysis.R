############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/CAIL/YuklProt3/", 
                  "~/gdrive/CAIL/YuklProt3/")
setwd(baseDir)

library(tidyverse)
library(ggrepel)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process text files ############
allPeptides <- read.delim("txt/allPeptides.txt", check.names = FALSE)
peptides <- read.delim("txt/peptides.txt", check.names = FALSE)
evidence <- read.delim("txt/evidence.txt", check.names = FALSE)
protGroups <- read.delim("txt/proteinGroups.txt", check.names = FALSE)
msmsScan <- read.delim("txt/msmsScans.txt", check.names = FALSE)
msms <- read.delim("txt/msms.txt", check.names = FALSE)
modPep <- read.delim("txt/modificationSpecificPeptides.txt", check.names = FALSE)

############ Protein processing, normalization, and DE ############
length(unique(protGroups[str_split(protGroups$`Protein IDs`, ";", simplify = TRUE)[,2] == "",]$`Protein IDs`)) # 1,381
allProts <- unique(c(as.matrix(str_split(protGroups$`Protein IDs`, ";", simplify = TRUE))))

# Removing those not present in at least 2 replicates and not present in 2 pooled:
protDF1 <- protGroups[, names(protGroups) %in% c("Protein IDs") | grepl("Intensity ", names(protGroups))]
pres <- apply(protDF1[, names(protDF1) != "Protein IDs"], 1, function(x) sum(x > 0)) > 1 
table(pres) # 1,285 pres vs 111
protDF1 <- protDF1[pres, ]
pres <- apply(protDF1[, grepl("pool", names(protDF1))], 1, function(x) sum(x > 0)) > 1
table(pres) # 1,085 pres vs 200
protDF1 <- protDF1[pres, ]

# SampMin imputation:
protSampMins <- apply(protDF1[, names(protDF1) != "Protein IDs"], 2, function(x) min(x[x > 0]))
for(j in 1:length(protSampMins)){
  protDF1[,names(protSampMins)[j]][protDF1[,names(protSampMins)[j]] == 0] <- protSampMins[j]
}

# VSN:
protDF2 <- protDF1[, names(protDF1) != "Protein IDs"]
rownames(protDF2) <- protDF1$`Protein IDs`
protDF3 <- vsn::justvsn(as.matrix(protDF2))
vsn::meanSdPlot(protDF3)
protDF3 <- as.data.frame(protDF3)
protDF3$`Protein IDs` <- rownames(protDF3)

# Some transformations:
protDF1L <- protDF1 %>% gather(key = "Sample", value = "intensity", -`Protein IDs`)
protDF1L$Sample <- gsub("Intensity ", "", protDF1L$Sample)
protDF1L$intensity <- log2(protDF1L$intensity)

protDF3L <- protDF3 %>% gather(key = "Sample", value = "intensity", -`Protein IDs`)
protDF3L$Sample <- gsub("Intensity ", "", protDF3L$Sample)

# Plot intensity distributions:
ggplot(protDF1L, aes(y = intensity, x = Sample)) + geom_boxplot() + ylim(0, 35)
ggplot(protDF3L, aes(y = intensity, x = Sample)) + geom_boxplot() + ylim(0, 35)

# Limma Fold Change differences:
protDF3$`Protein IDs` <- NULL
des1 <- model.matrix(~ 0 + factor(c(rep(1, 4), rep(2, 4))))
colnames(des1) <- c("hnox", "wt")
des1c <- limma::makeContrasts(hnox - wt, levels = des1)
protFit0 <- limma::lmFit(protDF3[, grepl("hnox|wt", names(protDF3))], design = des1)
protFit0 <- limma::eBayes(protFit0)
protFit1 <- limma::eBayes(limma::contrasts.fit(protFit0, des1c))
protDE <- data.frame(log2FC = protFit1$coefficients[,1], pValue = protFit1$p.value[,1])
protDE$protID <- rownames(protDE)
protDE$qValue <- qvalue::qvalue(protDE$pValue)$qvalues

# Volcano plot:
ggplot(protDE, aes(x = log2FC, y = -log10(pValue))) + geom_point() + 
  geom_hline(yintercept = -log10(.05), color = "red", lty = 2, lwd = 2)

# PCA:
protPCA <- prcomp(t(protDF3), center = TRUE, scale = TRUE)
protPCAx <- as.data.frame(protPCA$x[,1:3])
protPCAx$Sample <- gsub("Intensity ", "", rownames(protPCAx))
protPCAx$Phenotype <- gsub("\\d", "", protPCAx$Sample)
ggplot(protPCAx, aes(x = PC1, y = PC2, color = Phenotype, label = Sample)) + geom_point() + geom_text_repel()

# Agreement plots
ggplot(protDF3, aes(x = `Intensity hnox1`, y = `Intensity hnox2`)) + geom_point()
ggplot(protDF3, aes(x = `Intensity hnox3`, y = `Intensity hnox2`)) + geom_point()
ggplot(protDF3, aes(x = `Intensity hnox1`, y = `Intensity wt1`)) + geom_point()
cor(protDF3$`Intensity hnox1`, protDF3$`Intensity hnox2`)
cor(protDF3$`Intensity hnox3`, protDF3$`Intensity hnox2`)
cor(protDF3$`Intensity hnox1`, protDF3$`Intensity wt1`)
idk <- cor(protDF3)
diag(idk) <- .8
corrplot::corrplot(idk, is.corr = FALSE, diag = FALSE, type = "lower", method = "color", addCoef.col = "black")

############ Peptide processing and normalization ############
peptides2 <- as.data.frame(peptides[, names(peptides) == "id" | grepl("Intensity", names(peptides))])
peptides2$Intensity <- NULL

# Removing those not present in at least 2 replicates and not present in 2 pooled:
pres <- apply(peptides2[, names(peptides2) != "id"], 1, function(x) sum(x > 0)) > 1 
table(pres) # 5,162 pres vs 1,412
peptides2 <- peptides2[pres, ]
pres <- apply(peptides2[, grepl("pool", names(peptides2))], 1, function(x) sum(x > 0)) > 1
table(pres) # 3,588 pres vs 1,574
peptides2 <- peptides2[pres, ]

# SampMin imputation:
pepSampMins <- apply(peptides2[, names(peptides2) != "id"], 2, function(x) min(x[x > 0]))
for(j in 1:length(pepSampMins)){
  peptides2[,names(pepSampMins)[j]][peptides2[,names(pepSampMins)[j]] == 0] <- pepSampMins[j]
}

# Normalization and FC for present in both:
rownames(peptides2) <- peptides2$id
pepVsn <- vsn::justvsn(as.matrix(peptides2[, names(peptides2) != "id"]))
vsn::meanSdPlot(pepVsn)
peptides3 <- as.data.frame(pepVsn)

# PCA:
pepPCA <- prcomp(t(peptides3), center = TRUE, scale = TRUE)
pepPCAx <- as.data.frame(pepPCA$x[,1:3])
pepPCAx$Sample <- gsub("Intensity ", "", rownames(pepPCAx))
pepPCAx$Phenotype <- gsub("\\d", "", pepPCAx$Sample)
ggplot(pepPCAx, aes(x = PC1, y = PC2, color = Phenotype, label = Sample)) + geom_point() + geom_text_repel()

# Correlations:
idk <- cor(peptides3)
corrplot::corrplot(idk, is.corr = FALSE, diag = FALSE, type = "lower", method = "color", addCoef.col = "black")

############ Find a good peptide ############
peptides2$q20 <- apply(peptides2[, names(peptides2) != "id"], 1, function(x) quantile(x, probs = .20))
