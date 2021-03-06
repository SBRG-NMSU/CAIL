############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
library(tidyverse)
library(ggrepel)

oldPar <- par()
baseDir <- ifelse(Sys.info()["sysname"] == "Windows", "C:/Users/ptrainor/gdrive/CAIL/Energy/", "~/gdrive/CAIL/Energy/")
setwd(baseDir)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process ############
df1 <- readxl::read_excel("ht-modi-3.xlsx") %>% as.data.frame()
table(df1$`sample ID`)
df1$`Molecular Formula` <- gsub(" ", "", df1$`Molecular Formula`)

df1$dbeMF <- paste(df1$DBE, df1$`Molecular Formula`, sep = "_")
df2 <- df1 %>% group_by(dbeMF) %>% summarize(n = n()) %>% as.data.frame()
df3 <- df1 %>% group_by(`Molecular Formula`) %>% summarize(n = n())

# 1.     Ww12, ww14, ww17, ww15
# 2.     Ww15 vs 95-2
# 3.     Ww10 vs ww17
# 4.     86-14 vs ww17
# 5.     Ww10 vs 86-41
# 6.     86-46 vs ww17

############ Make wide datasets ############
# Presence vs absence wide:
df2$ww10 <- df2$ww12 <- df2$ww14 <- df2$ww15 <- df2$ww17ss1 <- df2$ww17ss2 <- df2$`86-41` <- df2$`86-46` <- df2$`95-2` <- 0
mapDF <- data.frame(sName = c("ww10", "ww12", "ww14", "ww15", "ww17ss1", "ww17ss2", "86-41", "86-46", "95-2"), 
                    lName = c("ww10-SA", "ww12-SS", "ww14 Bulk", "ww15 SA", "ww 17 SS1", "ww 17 SS 2", "HT62006-86-41", "HT62006-86-46", "HT62006-95-2"))
for(i in 1:nrow(mapDF)){
  temp1 <- df1[df1$`sample ID` == mapDF$lName[i],]
  df2[,mapDF$sName[i]] <- as.integer(df2$dbeMF %in% temp1$dbeMF)
}

# Intensity wide:
df3 <- df2[, names(df2) != "n"]
df3$ww10 <- df3$ww12 <- df3$ww14 <- df3$ww15 <- df3$ww17ss1 <- df3$ww17ss2 <- df3$`86-41` <- df3$`86-46` <- df3$`95-2` <- 0
for(i in 1:nrow(mapDF)){
  temp1 <- df1[df1$`sample ID` == mapDF$lName[i],]
  df3[,mapDF$sName[i]] <- temp1$`Rel. Abundance`[match(df3$dbeMF, temp1$dbeMF)]
}
rownames(df3) <- df3$dbeMF
df3$dbeMF <- NULL

# Super faint abundance:
df3 <- df3[apply(df3, 1, function(x) sum(x, na.rm = TRUE)) > 1,] # 18870 to 4802

# Present in one only:
df3 <- df3[apply(df3, 1, function(x) sum(!is.na(x)) > 2), ] # 4802 to 4138 

# Zeros:
df4 <- df3b <- df3
for(i in 1:nrow(df3)){
  for(j in 1:ncol(df3)){
    df3[i,j] <- ifelse(is.na(df3[i,j]), 0.00001, df3[i,j])
    df3b[i,j] <- ifelse(is.na(df3b[i,j]), 0.0, df3b[i,j])
    #df4[i,j] <- ifelse(is.na(df4[i,j]), 0.0, df4[i,j])
  }
}

############ Multivariate ############
# ZIFA: 
write.csv(t(df3b), file = "df3b.csv")
zifaScores <- read.csv("ZIFAScores.tab", sep = " ", header = FALSE)
zifaLoadings <- read.csv("ZIFALoad.tab", sep = " ", header = FALSE)
df3b <- cbind(df3b, zifaLoadings)
df3b$dbeMF <- rownames(df3b)

df3b$dbeMF <- rownames(df3b)
df3b$dbe <- as.integer(str_split(df3b$dbeMF, "_", simplify = TRUE)[,1])
df3b$MF <- str_split(df3b$dbeMF, "_", simplify = TRUE)[,2]
df3b$CNumb <- as.integer(gsub("C", "", str_extract(df3b$MF, "C\\d+")))

ggplot(df3b %>% filter(abs(V1) > 3), aes(y = dbe, x = CNumb, color = V1)) + geom_point() + 
  scale_color_gradient2() + 
  xlim(0, 40) + ylim(0, 25)

ggplot(df3b %>% filter(abs(V2) > 3), aes(y = dbe, x = CNumb, color = V2)) + geom_point() + 
  scale_color_gradient2() + 
  xlim(0, 40) + ylim(0, 25)

# Probabalistic PCA
m1 <- t(as.matrix(df4))
m2 <- scale(m1, center = FALSE, scale = FALSE)
pca2 <- pcaMethods::pca(m2, method = "ppca", nPcs = 4, seed = 3)
pca2L <- as.data.frame(pca2@loadings)
pca2X <- as.data.frame(pca2@scores)
pca2X$lab <- rownames(pca2X)
ggplot(pca2X, aes(x = PC1, y = PC2, label = lab)) + geom_point() + geom_label_repel()
ggplot(pca2X, aes(x = PC2, y = PC3, label = lab)) + geom_point() + geom_label_repel()

pca2L$dbeMF <- rownames(pca2L)
pca2L$dbe <- as.integer(str_split(pca2L$dbeMF, "_", simplify = TRUE)[,1])
pca2L$MF <- str_split(pca2L$dbeMF, "_", simplify = TRUE)[,2]
pca2L$CNumb <- as.integer(gsub("C", "", str_extract(pca2L$MF, "C\\d+")))

ggplot(pca2L %>% filter(abs(PC1) > .02), aes(y = dbe, x = CNumb, color = PC1 * 100)) + geom_point() + scale_color_distiller(palette = "Spectral") + 
  xlim(0, 40) + ylim(0, 25)
ggplot(pca2L %>% filter(abs(PC2) > .02), aes(y = dbe, x = CNumb, color = PC2 * 100)) + geom_point() + scale_color_distiller(palette = "Spectral") + 
  xlim(0, 40) + ylim(0, 25)
ggplot(pca2L %>% filter(abs(PC3) > .02), aes(y = dbe, x = CNumb, color = PC3 * 100)) + geom_point() + scale_color_distiller(palette = "Spectral") + 
  xlim(0, 40) + ylim(0, 25)

myPal <- colorRampPalette(c("red", "white", "blue"))(n = 299)
png(file = "./Plots/heatmap.png", height = 40, width = 10, units = "in", res = 800)
gplots::heatmap.2(log(as.matrix(df3)), cexRow = .08, trace = "none", density.info = "none", col = myPal)
dev.off()

png(file = "./Plots/heatmap2.png", height = 7, width = 7, units = "in", res = 500)
gplots::heatmap.2(log(as.matrix(df3)), labRow = "", trace = "none", density.info = "none", col = myPal)
dev.off()

sampOrder <- c("ww17ss2", "ww15", "ww14", "w10", "ww17ss1", "ww12", "86-41", "95-2", "86-46")

############ Weighted network analysis ############
df3$dbMFE <- rownames(df3)
df3w <- data.table::transpose(df3, keep.names = "samp", make.names = "dbMFE")
rownames(df3w) <- df3w$samp
df3w$samp <- NULL
df3w <- log2(df3w)
powers <- 1:20
sft <- WGCNA::pickSoftThreshold(df3w, powerVector = powers, verbose = 5, networkType = "signed hybrid", 
                                corOptions = list(use = "everything"))
par(mfrow = c(1, 2))
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.95, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = "red")
softPower <- 10
par(oldPar)

# Adjacency and distance matrices
adjacency <- WGCNA::adjacency(df3w, type = "signed hybrid", power = softPower, corOptions = list(use = "everything")) 
# Turn adjacency into topological overlap
TOM <- WGCNA::TOMsimilarity(adjacency, TOMType = "signed")
rownames(TOM) <- colnames(TOM) <- rownames(adjacency)
dissTOM <- 1 - TOM
tree <- hclust(as.dist(dissTOM), method = "complete")

# Module identification using dynamic tree cut:
dynamicMods <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = TRUE,
                                             minClusterSize = 20, method = "hybrid")
dynamicColors <- WGCNA::labels2colors(dynamicMods)

# Plot dendogram and module assignment:
png(filename = "Plots/WGCNA_Dendro1.png",height = 6, width = 10, units = "in", res = 600)
WGCNA::plotDendroAndColors(tree, cbind(dynamicColors), "Module Assignment", dendroLabels = FALSE, hang = 0.0, 
                           addGuide = TRUE, guideHang = 0.0, main = "Modules")
dev.off()

dynamicColors2 <- WGCNA::mergeCloseModules(df3w, dynamicColors, cutHeight = .01)$colors
WGCNA::plotDendroAndColors(tree, cbind(dynamicColors2), "Module Assignment", dendroLabels = FALSE, hang = 0.0, 
                           addGuide = TRUE, guideHang = 0.0, main = "Modules")

# Plot distnace matrix heatmap w/ module assignment
png(filename = "Plots/WGCNA_Heatmap1.png", height = 6, width = 6, units = "in", res = 600)
WGCNA::TOMplot(dissTOM, tree, dynamicColors, main = "Module heatmap")
dev.off()

# Module-feature mapping:
modDF <- data.frame(feature = rownames(TOM), module = dynamicColors)
modDF2 <- modDF
modDF2 <- modDF2 %>% left_join(df3b %>% select(-V1, -V2), by = c("feature" = "dbeMF"))
# writexl::write_xlsx(modDF2, path = "featuresWModules.xlsx")
as.data.frame(table(modDF$module)) %>% arrange(desc(Freq))

########### Module Plots ############
# turquoise module:
turquoiseTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "turquoise"], modDF$feature[modDF$module == "turquoise"]]),
                          method = "complete")
turquoiseExpression <- df3w[, modDF$feature[modDF$module == "turquoise"]]
png(filename = "Plots/modules/turquoise.png", height = 6, width = 24, units = "in", res = 600)
turquoiseColors <- WGCNA::numbers2colors(turquoiseExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = turquoiseTOMhClust, colors = t(turquoiseColors),
                           groupLabels = rownames(turquoiseExpression),
                           cex.dendroLabels = .17, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# blue module:
blueTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "blue"], modDF$feature[modDF$module == "blue"]]),
                          method = "complete")
blueExpression <- df3w[, modDF$feature[modDF$module == "blue"]]
png(filename = "Plots/modules/blue.png", height = 6, width = 20, units = "in", res = 600)
blueColors <- WGCNA::numbers2colors(blueExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = blueTOMhClust, colors = t(blueColors),
                           groupLabels = rownames(blueExpression),
                           cex.dendroLabels = .25, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# brown module:
brownTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "brown"], modDF$feature[modDF$module == "brown"]]),
                              method = "complete")
brownExpression <- df3w[, modDF$feature[modDF$module == "brown"]]
png(filename = "Plots/modules/brown.png", height = 6, width = 15, units = "in", res = 600)
brownColors <- WGCNA::numbers2colors(brownExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = brownTOMhClust, colors = t(brownColors),
                           groupLabels = rownames(brownExpression),
                           cex.dendroLabels = .3, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# yellow module:
yellowTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "yellow"], modDF$feature[modDF$module == "yellow"]]),
                         method = "complete")
yellowExpression <- df3w[, modDF$feature[modDF$module == "yellow"]]
png(filename = "Plots/modules/yellow.png", height = 6, width = 15, units = "in", res = 600)
yellowColors <- WGCNA::numbers2colors(yellowExpression, signed = TRUE, commonLim = FALSE)
WGCNA::plotDendroAndColors(dendro = yellowTOMhClust, colors = t(yellowColors),
                           groupLabels = rownames(yellowExpression),
                           cex.dendroLabels = .35, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# green module:
greenTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "green"], modDF$feature[modDF$module == "green"]]),
                        method = "average")
greenExpression <- df3w[, modDF$feature[modDF$module == "green"]]
png(filename = "Plots/modules/green.png", height = 6, width = 12, units = "in", res = 600)
greenColors <- WGCNA::numbers2colors(greenExpression, signed = TRUE, commonLim = TRUE)
WGCNA::plotDendroAndColors(dendro = greenTOMhClust, colors = t(greenColors),
                           groupLabels = rownames(greenExpression),
                           cex.dendroLabels = 0.5, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# red module:
redTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "red"], modDF$feature[modDF$module == "red"]]),
                         method = "average")
redExpression <- df3w[, modDF$feature[modDF$module == "red"]]
png(filename = "Plots/modules/red.png", height = 6, width = 12, units = "in", res = 600)
redColors <- WGCNA::numbers2colors(redExpression, signed = TRUE, commonLim = TRUE)
WGCNA::plotDendroAndColors(dendro = redTOMhClust, colors = t(redColors),
                           groupLabels = rownames(redExpression),
                           cex.dendroLabels = 0.5, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# black module:
blackTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "black"], modDF$feature[modDF$module == "black"]]),
                       method = "average")
blackExpression <- df3w[, modDF$feature[modDF$module == "black"]]
png(filename = "Plots/modules/black.png", height = 6, width = 10, units = "in", res = 600)
blackColors <- WGCNA::numbers2colors(blackExpression, signed = TRUE, commonLim = TRUE)
WGCNA::plotDendroAndColors(dendro = blackTOMhClust, colors = t(blackColors),
                           groupLabels = rownames(blackExpression),
                           cex.dendroLabels = 0.5, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# magenta module:
magentaTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "magenta"], modDF$feature[modDF$module == "magenta"]]),
                         method = "average")
magentaExpression <- df3w[, modDF$feature[modDF$module == "magenta"]]
png(filename = "Plots/modules/magenta.png", height = 6, width = 10, units = "in", res = 600)
magentaColors <- WGCNA::numbers2colors(magentaExpression, signed = TRUE, commonLim = TRUE)
WGCNA::plotDendroAndColors(dendro = magentaTOMhClust, colors = t(magentaColors),
                           groupLabels = rownames(magentaExpression),
                           cex.dendroLabels = 0.5, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

# pink module:
pinkTOMhClust <- hclust(as.dist(TOM[modDF$feature[modDF$module == "pink"], modDF$feature[modDF$module == "pink"]]),
                           method = "average")
pinkExpression <- df3w[, modDF$feature[modDF$module == "pink"]]
png(filename = "Plots/modules/pink.png", height = 6, width = 10, units = "in", res = 600)
pinkColors <- WGCNA::numbers2colors(pinkExpression, signed = TRUE, commonLim = TRUE)
WGCNA::plotDendroAndColors(dendro = pinkTOMhClust, colors = t(pinkColors),
                           groupLabels = rownames(pinkExpression),
                           cex.dendroLabels = 0.5, cex.colorLabels = .7, marAll = c(1, 6, 3, 1), main = "")
dev.off()

########### Module PCs ###########
mEigen1 <- WGCNA::moduleEigengenes(df3w, dynamicColors, impute = FALSE, nPC = 1, align = "along average", 
                                   excludeGrey = TRUE, grey = if (is.numeric(colors)) 0 else "grey", 
                                   softPower = 10, scale = TRUE, verbose = 5, indent = 1)
mEigen1T <- as.data.frame(t(mEigen1$eigengenes))
mEigen1T$module <- gsub("ME", "", rownames(mEigen1T))
writexl::write_xlsx(mEigen1T, path = "Eigen.xlsx")
