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
man1 <- readxl::read_excel("Swine-Manure-MF.xlsx", sheet = 1)
man2 <- readxl::read_excel("Swine-Manure-MF.xlsx", sheet = 2)
man3 <- readxl::read_excel("Swine-Manure-MF.xlsx", sheet = 3)
man <- rbind(man1, man2, man3)
rm(man1, man2, man3)

sludge1 <- readxl::read_excel("CCCSD-Sludge-MF.xlsx", sheet = 1)
sludge2 <- readxl::read_excel("CCCSD-Sludge-MF.xlsx", sheet = 2)
sludge3 <- readxl::read_excel("CCCSD-Sludge-MF.xlsx", sheet = 3)
sludge <- rbind(sludge1, sludge2, sludge3)
rm(sludge1, sludge2, sludge3)

df1 <- rbind(man, sludge)

# Get rid of duplicates:
df1$uniqueCheck <- paste(df1$`Samle ID`, df1$MF, df1$DBE, df1$RA, sep = "_")
df1$uniqueCheck2 <- duplicated(df1$uniqueCheck)
table(df1$`Samle ID`, df1$uniqueCheck2)
df1 <- df1[!df1$uniqueCheck2,] # 21,975 to 21,970 # double check this
df1 <- df1[, !grepl("unique", names(df1))]
df1$mFDBE <- paste(df1$MF, df1$DBE, sep = "_")

# Long to wide for imputation:
annoDF <- df1 %>% select(mFDBE, Class, MF, DBE) %>% distinct()
df1 <- df1 %>% select(SampleID = `Samle ID`, RA, mFDBE)

# Temporarily remove non-unique:
# weird ones:
df1[c(1934, 1935, 5025, 8667, 5211, 8853, 6210, 9852),]
df1 <- df1[-c(1934, 1935, 5025, 8667, 5211, 8853, 6210, 9852),]
df1w <- df1 %>% spread(key = mFDBE, value = RA, fill = 0)
df1 <- df1w %>% gather(key = "mFDBE", value = RA, -SampleID)

############ Filtering, Imputation & Normalization ############
# Filter those that aren't present in 2 or more samples:
df1 <- df1 %>% left_join(df1 %>% group_by(mFDBE) %>% summarize(sumNonZeroInd = sum(RA > 0) > 1))
length(unique(df1$mFDBE)) # 3,922
df1 <- df1 %>% filter(sumNonZeroInd)
length(unique(df1$mFDBE)) # 3,901

# Minimum / 2 imputation:
minNonZero <- function(x){
  x <- x[x > 0]
  return(min(x))
}
df1 <- df1 %>% left_join(df1 %>% group_by(mFDBE) %>% summarize(mNZ = minNonZero(RA)))
df1$RA[df1$RA == 0] <- df1$mNZ[df1$RA == 0]
df1 <- df1 %>% select(-sumNonZeroInd, -mNZ)

# Log transform:
df1$Abundance <- log2(df1$RA)

# Add group:
df1$Type <- gsub(" \\d", "", df1$SampleID)
table(df1$Type)

# Histogram of abundances:
ggplot(df1, aes(x = SampleID, y = Abundance)) + geom_boxplot()

# Is normalization needed?

############ Statistical Test ############
uniqueFeatures <- unique(df1$mFDBE)
df2 <- data.frame(feature = uniqueFeatures, logFC = NA, pValue = NA)
for(i in 1:nrow(df2)){
  temp1 <- df1 %>% filter(mFDBE == df2$feature[i])
  temp2 <- t.test(Abundance ~ Type, data = temp1)
  df2$logFC[i] <- temp2$estimate[1] - temp2$estimate[2]
  df2$pValue[i] <- temp2$p.value
}
df2$qValue <- qvalue::qvalue(df2$pValue)$qvalue

ggplot(df2, aes(x = logFC, y = -log10(pValue))) + geom_point()

############ Heatmap ############
df1w <- as.data.frame(df1w)
rownames(df1w) <- df1w$SampleID
df1w$SampleID <- NULL
heatmap(as.matrix(df1w))

plot(as.numeric(df1w[1,]), as.numeric(df1w[2,]))
plot(as.numeric(df1w[2,]), as.numeric(df1w[3,]), ylim = c(0, 25))
plot(as.numeric(df1w[1,]), as.numeric(df1w[3,]), ylim = c(0, 25))
cor(t(df1w))
corrplot::corrplot(cor(t(df1w)))

plot(as.numeric(df1w[4,]), as.numeric(df1w[5,]))
plot(as.numeric(df1w[5,]), as.numeric(df1w[6,]))
plot(as.numeric(df1w[4,]), as.numeric(df1w[6,]), ylim = c(0, 14))
