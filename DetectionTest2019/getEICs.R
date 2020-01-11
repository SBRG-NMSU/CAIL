options(stringsAsFactors = FALSE)

library(tidyverse)
library(xcms)

setwd("~/CAIL/DetectionTest")

files <- list.files(pattern = ".mzML")
files <- files[(grepl("Blank", files) | grepl("ppb", files) ) & grepl("MS1", files)]
df0 <- readMSData(files = files, msLevel. = 1)

# register(bpstart(SnowParam(2)))
register(SerialParam())

bpis <- chromatogram(df0, aggregationFun = "max")
tics <- chromatogram(df0, aggregationFun = "sum")

ticsList <- lapply(tics, function(x) data.frame(rtime = x@rtime, intensity = x@intensity, name = names(x@rtime)[1]))
ticsDF <- do.call("rbind", ticsList)

nameDF <- data.frame(old = unique(ticsDF$name), Name = gsub(".mzML", "", files))
ticsDF <- ticsDF %>% left_join(nameDF, by = c("name" = "old"))

png(file = "Images/PPBTics.png", height = 4, width = 8, units = "in", res = 600)
ggplot(ticsDF, aes(x = rtime, y = intensity, group = Name, color = Name)) + geom_line() + theme_bw() + 
  labs(y = "TIC", x = "RT")
dev.off()

standards <- readxl::read_excel('InternalStandards.xlsx') %>% as.data.frame()
standards$da <- standards$`M+H`
standards$lower <- standards$da - standards$da * 5 / 1000000
standards$upper <- standards$da + standards$da * 5 / 1000000

allPeaks <- data.frame()
for(i in 1:nrow(standards)){
  eic1 <- chromatogram(df0, mz = c(standards$lower[i], standards$upper[i]), missing = 100)
  peaks <- findChromPeaks(eic1, CentWaveParam(ppm = 10, snthresh = 10, peakwidth = c(15, 70), prefilter = c(1, 5000)))
  peaks2 <- lapply(peaks, function(x) x@chromPeaks)
  names(peaks2) <- peaks@phenoData@data$sampleNames
  for(j in 1:length(peaks2)){
    peaks2[[j]] <- as.data.frame(peaks2[[j]])
    if(nrow(peaks2[[j]]) == 0) peaks2[[j]] <- peaks2[[j]][1,]
    peaks2[[j]]$fName <- names(peaks2)[j]
  }
  peaks3 <- do.call("rbind", peaks2)
  peaks3$compound <- standards$Compound[i]
  png(file = paste0("Images/iSTD/PPB/iSTD_", standards$Compound[i], ".png"), height = 5, width = 7, res = 300, units = "in")
  plot(peaks, col = c("red", "green", "purple", "lightblue"), 
       main = paste0(standards$Compound[i], " ", round(standards$lower[i], 5), " ", round(standards$upper[i], 5)))
  dev.off()
  
  allPeaks <- rbind(allPeaks, peaks3)
}


