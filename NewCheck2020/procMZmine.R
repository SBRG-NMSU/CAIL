############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 600)
oldPar <- par()

setwd("C:/Users/ptrainor/gdrive/CAIL/NewCheck2020")

# Import mzmine2 data:
small <- read.csv("NanomateSmallQuan.csv")
names(small) <- gsub("_", "-", names(small))

peakData <- small[ , grepl("raw.Peak.FWHM|.area|.m.z|tailing.factor|asymmetry.factor|.RT", names(small)) &
                     !grepl("row", names(small))]
names(peakData) <- gsub("raw.Peak", "", names(peakData))
names(peakData) <- gsub("\\.", "_", names(peakData))

anno <- small[ , grepl("row.ID|row.m.z|row.retention.time|row.identity..main.ID", names(small)) & 
                 !grepl("details", names(small))]
names(anno)[match(c("row.ID","row.identity..main.ID.", "row.retention.time", "row.m.z"), names(anno))] <- 
  c("CmpdID", "CmpdName", "RTavg", "MZavg")
anno <- anno[, c("CmpdID", "CmpdName", "RTavg", "MZavg")]

peakData$CmpdID <- anno$CmpdID

missingP <- apply(peakData[, grepl("m_z_max", names(peakData))], 1, function(x) sum(ifelse(x == 0, 1, 0))) / nrow(peakData)
anno$missingP <- missingP

x <- peakData[, grepl("area", names(peakData))]
rsdFun <- function(x){
  x[x == 0] <- NA
}
vars <- c("m_z_min")
peakData[, grepl(vars[1], names(peakData))]
