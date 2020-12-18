############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/CAIL/YuklProt2/", 
                  "~/gdrive/CAIL/YuklProt2/")
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
length(unique(protGroups[str_split(protGroups$`Protein IDs`, ";", simplify = TRUE)[,2] == "",]$`Protein IDs`)) # 1,655
allProts <- unique(c(as.matrix(str_split(protGroups$`Protein IDs`, ";", simplify = TRUE))))

# Removing those not present in at least 2 replicates:
protDF1 <- protGroups[, names(protGroups) %in% c("Protein IDs") | grepl("Intensity ", names(protGroups))]
pres <- apply(protDF1[, names(protDF1) != "Protein IDs"], 1, function(x) sum(x > 0)) > 1 # 1,471
protDF1 <- protDF1[pres, ]

# SampMin imputation:
protSampMins <- apply(protDF1[, names(protDF1) != "Protein IDs"], 2, function(x) min(x[x > 0]))
for(j in 2:ncol(protDF1)){
  protDF1[,j][protDF1[,j] == 0] <- protSampMins[match(names(protDF1)[j], names(protSampMins))]
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
des1 <- model.matrix(~ 0 + factor(c(rep(1, 2), rep(2, 2))))
colnames(des1) <- c("hnox", "wt")
des1c <- limma::makeContrasts(hnox - wt, levels = des1)
protFit0 <- limma::lmFit(protDF3[,c(2,3,4,6)], design = des1)
protFit0 <- limma::eBayes(protFit0)
protFit1 <- limma::eBayes(limma::contrasts.fit(protFit0, des1c))
protDE <- data.frame(log2FC = protFit1$coefficients[,1], pValue = protFit1$p.value[,1])
protDE$protID <- rownames(protDE)
protDE$qValue <- qvalue::qvalue(protDE$pValue)$qvalues

# Volcano plot:
ggplot(protDE, aes(x = log2FC, y = -log10(qValue))) + geom_point() + 
  geom_hline(yintercept = -log10(.05), color = "red", lty = 2, lwd = 2)

# PCA:
protPCA <- prcomp(t(protDF3), center = TRUE, scale = TRUE)
protPCAx <- as.data.frame(protPCA$x[,1:3])
protPCAx$Sample <- gsub("Intensity ", "", rownames(protPCAx))
protPCAx$Phenotype <- gsub("\\d", "", protPCAx$Sample)
ggplot(protPCAx, aes(x = PC1, y = PC2, color = Phenotype, label = Sample)) + geom_point() + geom_text_repel()

# Sensitivty analysis
protPCA2 <- prcomp(t(protDF3[,-1]), center = TRUE, scale = TRUE)
protPCAx2 <- as.data.frame(protPCA2$x[,1:3])
protPCAx2$Sample <- gsub("Intensity ", "", rownames(protPCAx2))
protPCAx2$Phenotype <- gsub("\\d", "", protPCAx2$Sample)
ggplot(protPCAx2, aes(x = PC1, y = PC2, color = Phenotype, label = Sample)) + geom_point() + geom_text_repel()

ggplot(protDF3, aes(x = `Intensity hnox1`, y = `Intensity hnox2`)) + geom_point()
ggplot(protDF3, aes(x = `Intensity hnox3`, y = `Intensity hnox2`)) + geom_point()
ggplot(protDF3, aes(x = `Intensity hnox1`, y = `Intensity wt1`)) + geom_point()
cor(protDF3$`Intensity hnox1`, protDF3$`Intensity hnox2`)
cor(protDF3$`Intensity hnox3`, protDF3$`Intensity hnox2`)
cor(protDF3$`Intensity hnox1`, protDF3$`Intensity wt1`)
idk <- cor(protDF3)
diag(idk) <- .8
corrplot::corrplot(idk, cl.lim = c(0.75, .87), is.corr = FALSE, diag = FALSE, type = "lower")

# Get peptides for those:
pG3pep <- peptides[peptides$id %in% do.call("c", sapply(protGroups[protGroups$id %in% protGroups3$id,]$`Peptide IDs`, function(x) str_split(x, ";"))),]
pG3pepE <- evidence[evidence$`Peptide ID` %in% do.call("c", sapply(protGroups[protGroups$id %in% protGroups3$id,]$`Peptide IDs`, function(x) str_split(x, ";"))),]

# Protein groups vs proteins:
length(unique(protGroups[str_split(protGroups$`Protein IDs`, ";", simplify = TRUE)[,2] == "",]$`Protein IDs`))

# All protein DE:
protGroups2b <- protGroups %>% select(id, `Intensity WT`, `Intensity HNOX`, `Protein IDs`,
                                      `Fasta headers`)
protVsnb <- vsn::justvsn(as.matrix(protGroups2b[, 2:3]))
protGroups2b$intWT <- protVsnb[,"Intensity WT"]
protGroups2b$intHNOX <- protVsnb[,"Intensity HNOX"]
protGroups2b$FC <- protGroups2b$intHNOX - protGroups2b$intWT
writexl::write_xlsx(protGroups2b, path = "Contrasts_20200910.xlsx")

############ Peptide processing and normalization ############
# Present vs. absent:
peptides2 <- peptides %>% select(id, `Intensity WT`, `Intensity HNOX`) %>% as.data.frame()
peptides2$presWT <- ifelse(peptides2$`Intensity WT` > 0, TRUE, FALSE)
peptides2$presHNOX <- ifelse(peptides2$`Intensity HNOX` > 0, TRUE, FALSE)
peptidesPresAbs <- xtabs(~ presHNOX + presWT, data = peptides2)

# Normalization and FC for present in both:
peptides2 <- peptides2 %>% filter(presWT & presHNOX)
pepVsn <- vsn::justvsn(as.matrix(peptides2[, 2:3]))
peptides2$IntWT <- pepVsn[,"Intensity WT"]
peptides2$IntHNOX <- pepVsn[,"Intensity HNOX"]
peptides2$FC <- peptides2$IntHNOX - peptides2$IntWT
peptides2$FCType <- ifelse(peptides2$FC > 1, "Higher", ifelse(peptides2$FC < -1, "Lower", "Ind"))
peptidesFC <- xtabs(~ FCType, data = peptides2)
