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
df3 <- df3[apply(df3, 1, function(x) sum(x, na.rm = TRUE)) > 1,]

# Present in one only:
df3 <- df3[apply(df3, 1, function(x) sum(!is.na(x)) > 2), ]

# Zeros:
df4 <- df3b <- df3
for(i in 1:nrow(df3)){
  for(j in 1:ncol(df3)){
    df3[i,j] <- ifelse(is.na(df3[i,j]), 0.01, df3[i,j])
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

ggplot

# Regular PCA with small imputation:
pca1 <- prcomp(t(as.matrix(df3)), center = TRUE, scale = FALSE)
pca1X <- as.data.frame(pca1$x[,1:4])
pca1L <- as.data.frame(pca1$rotation[,1:4])
pca1X$lab <- rownames(pca1X)
ggplot(pca1X, aes(x = PC1, y = PC2, label = lab)) + geom_point() + geom_label_repel()
ggplot(pca1X, aes(x = PC2, y = PC3, label = lab)) + geom_point() + geom_label_repel()

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
png(file = "heatmap.png", height = 40, width = 10, units = "in", res = 800)
gplots::heatmap.2(log(as.matrix(df3)), cexRow = .08, trace = "none", density.info = "none", col = myPal)
dev.off()

png(file = "heatmap2.png", height = 7, width = 7, units = "in", res = 500)
gplots::heatmap.2(log(as.matrix(df3)), labRow = "", trace = "none", density.info = "none", col = myPal)
dev.off()

bdist <- vegan::vegdist(t(as.matrix(df3)), method = "bray")
fit1 <- MASS::isoMDS(bdist, k = 2)
plot(fit1$points)
