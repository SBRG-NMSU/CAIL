############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/CAIL/YuklProt4/", 
                  "~/gdrive/CAIL/YuklProt4/")
setwd(baseDir)

library(tidyverse)
library(ggrepel)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process text files ############
# allPeptides <- read.delim("txt/allPeptides.txt", check.names = FALSE)
# allPeptides$id_allPeptides <- 1:nrow(allPeptides)
# peptides <- read.delim("txt/peptides.txt", check.names = FALSE)
# evidence <- read.delim("txt/evidence.txt", check.names = FALSE)
# protGroups <- read.delim("txt/proteinGroups.txt", check.names = FALSE)
# msmsScan <- read.delim("txt/msmsScans.txt", check.names = FALSE)
# msms <- read.delim("txt/msms.txt", check.names = FALSE)
# modPep <- read.delim("txt/modificationSpecificPeptides.txt", check.names = FALSE)
# matchedFeatures <- read.delim("txt/matchedFeatures.txt", check.names = FALSE)
# depPeptides <- read.delim("txt/dependentPeptides.txt", check.names = FALSE)
# save.image("txt/txtData.RData")
load("txt/txtData.RData")

############ Protein processing, normalization, and DE ############
length(unique(protGroups[str_split(protGroups$`Protein IDs`, ";", simplify = TRUE)[,2] == "",]$`Protein IDs`)) # 1,075
allProts <- unique(c(as.matrix(str_split(protGroups$`Protein IDs`, ";", simplify = TRUE))))

## Filtering
# Removing those not present in at least 2 replicates and not present in 1 pooled:
protDF1 <- protGroups[, names(protGroups) %in% c("Protein IDs") | grepl("Intensity ", names(protGroups))]
pres <- apply(protDF1[, names(protDF1) != "Protein IDs"], 1, function(x) sum(x > 0)) > 1 
table(pres) # 1070 vs 14
protDF1 <- protDF1[pres, ]
pres <- apply(protDF1[, grepl("pool", names(protDF1))], 1, function(x) sum(x > 0)) > 0
table(pres) # 982 pres vs 88 # Should follow up on some of these
protDF1 <- protDF1[pres, ]

## Filtering con and reverse:
protDF1 <- protDF1 %>% filter(!`Protein IDs` %in% protDF1$`Protein IDs`[grepl("CON_", protDF1$`Protein IDs`)]) # 982 -> 977
protDF1 <- protDF1 %>% filter(!`Protein IDs` %in% protDF1$`Protein IDs`[grepl("REV_", protDF1$`Protein IDs`)]) # 977 -> 968

# SampMin imputation:
protSampMins <- apply(protDF1[, names(protDF1) != "Protein IDs"], 1, function(x) min(x[x > 0]))
protSampMinsDiv2 <- protSampMins / 2
for(i in 1:nrow(protDF1)){
  protDF1[i, names(protDF1) != "Protein IDs"][protDF1[i, names(protDF1) != "Protein IDs"] == 0] <- protSampMinsDiv2[i]
}

# VSN: # Why use VSN? https://academic.oup.com/bib/article/19/1/1/2562889
protDF2 <- protDF1[, names(protDF1) != "Protein IDs"]
rownames(protDF2) <- protDF1$`Protein IDs`
protDF3 <- vsn::justvsn(as.matrix(protDF2))
png(filename = "Plots/vsn.png", height = .75*6, width = .75*7, units = "in", res = 600)
vsn::meanSdPlot(protDF3)
dev.off()
protDF3 <- as.data.frame(protDF3)
protDF3$`Protein IDs` <- rownames(protDF3)
names(protDF1) <- gsub("wt", "WT", gsub("pool", "Pool", gsub("hnox", "H-NOX", gsub("Intensity ", "", names(protDF1)))))
names(protDF3) <- gsub("wt", "WT", gsub("pool", "Pool", gsub("hnox", "H-NOX", gsub("Intensity ", "", names(protDF3)))))

# Some transformations:
protDF1L <- protDF1 %>% gather(key = "Sample", value = "intensity", -`Protein IDs`)
protDF1L$Sample <- gsub("Intensity ", "", protDF1L$Sample)
protDF1L$intensity <- log2(protDF1L$intensity)

protDF3L <- protDF3 %>% gather(key = "Sample", value = "intensity", -`Protein IDs`)
protDF3L$Sample <- gsub("Intensity ", "", protDF3L$Sample)

# Plot intensity distributions:
p1 <- ggplot(protDF1L, aes(y = intensity, x = Sample)) + geom_boxplot() + ylim(0, 35) + ylab("Intensity (Log2 Scale)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggtitle("Pre-Normalization")
p2 <- ggplot(protDF3L, aes(y = intensity, x = Sample)) + geom_boxplot() + ylim(0, 35) + ylab("Normalized Intensity (Log2 Scale)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + ggtitle("Post-Normalization")
png("Plots/Norm.png", height = 4, width = 7, units = "in", res = 300)
show(gridExtra::grid.arrange(p1, p2, ncol = 2))
dev.off()

# Limma Fold Change differences:
protDF3$`Protein IDs` <- NULL
des1 <- model.matrix(~ 0 + factor(c(rep(1, 4), rep(2, 4))))
colnames(des1) <- c("hnox", "wt")
des1c <- limma::makeContrasts(hnox - wt, levels = des1)
protFit0 <- limma::lmFit(protDF3[, grepl("H-NOX|WT", names(protDF3))], design = des1)
protFit0 <- limma::eBayes(protFit0)
protFit1 <- limma::eBayes(limma::contrasts.fit(protFit0, des1c))
protDE <- data.frame(log2FC = protFit1$coefficients[,1], pValue = protFit1$p.value[,1])
protDE$protID <- rownames(protDE)
protDE$qValue <- qvalue::qvalue(protDE$pValue)$qvalues

# Volcano plot:
protDE$sigAnno <- ""
protDE$sigAnno[protDE$pValue < 0.01 & abs(protDE$log2FC) > 1] <- protDE$protID[protDE$pValue < 0.01 & abs(protDE$log2FC) > 1]
png(filename = "Plots/Volcano.png", height = 4, width = 5, units = "in", res = 2100)
set.seed(33)
ggplot(protDE, aes(x = log2FC, y = -log10(pValue))) + 
  geom_point(fill = "dodgerblue", color = "grey20", pch = 21, size = 1, alpha = .75) + 
  geom_text_repel(aes(x = log2FC, y = -log10(pValue), label = sigAnno), size = 2, box.padding = .5) + 
  geom_hline(yintercept = -log10(.01), color = "indianred", lty = 2, lwd = 1, alpha = .75) +
  geom_vline(xintercept = 1, color = "darkblue", lty = 2, lwd = 1, alpha = .75) +
  geom_vline(xintercept = -1, color = "darkblue", lty = 2, lwd = 1, alpha = .75) +
  xlim(-5, 5) + ylab("-Log10(p-value)") + xlab("Log2(Fold Change)")
dev.off()

# PCA:
protPCA <- prcomp(t(protDF3), center = TRUE, scale = TRUE)
protPCAx <- as.data.frame(protPCA$x[,1:3])
protPCAx$Sample <- rownames(protPCAx)
protPCAx$Phenotype <- gsub("\\d", "", protPCAx$Sample)
png(filename = "Plots/PCA.png", height = 4, width = 5, units = "in", res = 600)
ggplot(protPCAx, aes(x = PC1, y = PC2, color = Phenotype, label = Sample)) + geom_point() + geom_text_repel()
dev.off()
summary(protPCA)

# Agreement plots
png(filename = "Plots/corExample.png", height = 5, width = 5, units = "in", res = 300)
ggplot(protDF3, aes(x = `H-NOX1`, y = `H-NOX2`)) + geom_point(fill = "dodgerblue", color = "grey20", pch = 21, size = 1, alpha = .75)
dev.off()
ggplot(protDF3, aes(x = `H-NOX3`, y = `H-NOX2`)) + geom_point()
ggplot(protDF3, aes(x = `H-NOX1`, y = `WT1`)) + geom_point()
cor(protDF3$`H-NOX1`, protDF3$`H-NOX2`)
cor(protDF3$`H-NOX3`, protDF3$`H-NOX2`)
cor(protDF3$`H-NOX1`, protDF3$`WT1`)
corMat <- cor(protDF3)
diag(corMat) <- 1
rownames(corMat) <- colnames(corMat) <- gsub("wt", "WT", gsub("pool", "Pool", gsub("hnox", "H-NOX", gsub("Intensity ", "", rownames(corMat)))))
png(filename = "Plots/corMat.png", height = 5, width = 5, units = "in", res = 300)
corrplot::corrplot(corMat, is.corr = FALSE, diag = TRUE, method = "color", addCoef.col = "black", number.cex = .65)
dev.off()

# Heatmap:
png(filename = "Plots/heatmap_noLab.png", height = 8, width = 6, units = "in", res = 300)
pheatmap::pheatmap(as.matrix(protDF3), scale = "row", show_rownames = FALSE)
dev.off()
png(filename = "Plots/heatmap2.png", height = 10, width = 5, units = "in", res = 300)
pheatmap::pheatmap(as.matrix(protDF3), scale = "row", fontsize_row = 3)
dev.off()
heatmap(t(as.matrix(protDF3)))

# Join annotation data with DE data:
protDE$`Protein IDs` <- rownames(protDE)
protDE <- protDE %>% left_join(
  protGroups %>% select(`Protein IDs`:`Number of proteins`, `Sequence coverage [%]`, `Sequence length`, `Q-value`))

# Add a nice name:
protDE$ProtName <- unite(as.data.frame(str_split(str_split(protDE$`Fasta headers`, " OS=", simplify = TRUE)[,1],
                                                    "\\|", simplify = TRUE)[,c(1,3)]), "ProtName")[,1]
protDE$ProtName <- gsub("_PARDP", "", protDE$ProtName)

# Export data.frame:
writexl::write_xlsx(protDE,"Results/ProtDE.xlsx")

############ Peptides for some top proteins ############
protDETop <- protDE %>% filter(qValue < .10)
protDETop <- protDETop %>% arrange(pValue)
pepDETop1 <- peptides %>% filter(Proteins == protDETop$protID[1])
pepDETop1_1 <- peptides %>% filter(Sequence == "ISVLIDRK")
pepDETop1_2 <- peptides %>% filter(Sequence == "YMELAVQALHDAFELEKP")
evidence1_1 <- evidence %>% filter(id %in% c(29418, 29419, 29420, 29421, 29422, 29423, 29424, 29425))
evidence1_2 <- evidence %>% filter(id %in% c(66709:66725))

############ Plots ############
for(i in 1:nrow(protDETop)){
  temp1 <- protDETop$protID[i]
  pDETop1 <- data.frame(t(protDF3[temp1,]))
  names(pDETop1) <- "Intensity"
  pDETop1$Phenotype <- rownames(pDETop1)
  pDETop1$Phenotype <- gsub("Intensity ", "", pDETop1$Phenotype)
  pDETop1$Replicate <- parse_number(pDETop1$Phenotype)
  pDETop1$Phenotype <- gsub("\\d", "", pDETop1$Phenotype)
  pDETop1$Phenotype <- gsub("wt", "WT", gsub("hnox", "H-NOX", pDETop1$Phenotype))
  stat1 <- paste0("Log2 FC: ", round(protDETop$log2FC[i],3), "; Q-Value = ", round(protDETop$qValue[i], 3))
  png(file = paste0("Plots/DEProteins/", "DEProtein_", protDETop$protID[i], ".png"), height = 4, width = 4, units = "in", res = 300)
  show(ggplot(pDETop1 %>% filter(Phenotype != "pool"), aes(x = Phenotype, y = Intensity)) + 
    stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", color="red", width = 0.1) +
    stat_summary(fun.y = mean, geom = "point", color = "red", size = 3) + geom_point() + 
    ylab("Normalized Intensity (Log2 scale)") + ggtitle(protDETop$ProtName[i], stat1) +
    theme(plot.title = element_text(hjust = 0, size = 10), plot.subtitle = element_text(hjust = 0, size = 8)))
  dev.off()
}

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
