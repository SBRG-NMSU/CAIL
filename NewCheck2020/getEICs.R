options(stringsAsFactors = FALSE, scipen = 900)

library(tidyverse)
library(xcms)

setwd("~/CAIL/NewData/Nanomate_SmallColumn")
files <- list.files(pattern = ".mzML")
#files <- files[4:5]
ptm <- proc.time()
df0 <- readMSData(files = files, mode = "onDisk")
proc.time() - ptm
setwd("C:/Users/ptrainor/gdrive/CAIL/NewCheck")

#register(bpstart(SnowParam(4)))
register(SerialParam())

#bpis <- chromatogram(df0, aggregationFun = "max")

ptm <- proc.time()
tics <- chromatogram(df0, aggregationFun = "sum")
proc.time() - ptm

ticsList <- lapply(tics, function(x) data.frame(rtime = x@rtime, intensity = x@intensity, name = names(x@rtime)[1]))
ticsDF <- do.call("rbind", ticsList)

nameDF <- data.frame(old = unique(ticsDF$name), 
                     Name = gsub("nanomate_H2O_IS50ppb_PharMix1\\+2-5ppb_PesticideMix1\\+2-5ppb_SPE_Rec_20in_2ul_0p1Ijset-","",
                                 gsub(".mzML", "", files)), oldOld = files)
ticsDF <- ticsDF %>% left_join(nameDF, by = c("name" = "old"))

png(file = "Images/TicsAll.png", height = 2*4, width = 2*8, units = "in", res = 300)
ggplot(ticsDF, aes(x = rtime, y = intensity, group = Name, color = Name)) + geom_line(lwd = .25) + theme_bw() + 
  labs(y = "TIC", x = "RT")
dev.off()

png(file = "Images/TicsEveryThird.png", height = 2*4, width = 2*8, units = "in", res = 300)
ggplot(ticsDF %>% filter(Name %in% nameDF$Name[seq(from = 1, to = 15, by = 3)]), aes(x = rtime, y = intensity, group = Name, color = Name)) + 
  geom_line() + theme_bw() + 
  labs(y = "TIC", x = "RT")
dev.off()

standards <- readxl::read_excel('CompoundList.xlsx') %>% as.data.frame()
standards$da <- standards$MpH
standards$lower <- standards$da - standards$da * 3 / 1000000
standards$upper <- standards$da + standards$da * 3 / 1000000

allPeaks <- data.frame()
for(i in 1:nrow(standards)){
  print(i)
  eic1 <- chromatogram(filterFile(df0, seq(from = 1, to = 15, by = 3)), mz = c(standards$lower[i], standards$upper[i]), missing = 100)
  peaks <- findChromPeaks(eic1, CentWaveParam(ppm = 10, snthresh = 10, peakwidth = c(15, 40), prefilter = c(1, 5000)))
  peaks2 <- lapply(peaks, function(x) x@chromPeaks)
  names(peaks2) <- peaks@phenoData@data$sampleNames
  for(j in 1:length(peaks2)){
    peaks2[[j]] <- as.data.frame(peaks2[[j]])
    if(nrow(peaks2[[j]]) == 0) peaks2[[j]] <- peaks2[[j]][1,]
    peaks2[[j]]$fName <- names(peaks2)[j]
  }
  peaks3 <- do.call("rbind", peaks2)
  peaks3$compound <- standards$Compound[i]
  peaks3 <- peaks3 %>% left_join(nameDF, by = c("fName" = "oldOld"))
  png(file = paste0("Images/iSTD_", standards$Compound[i], "2.png"), height = 5, width = 7, res = 300, units = "in")
  plot(filterRt(peaks[1,1:3], rt = c(750, 850)), col = c("green", "red", "blue"),
       main = paste0(standards$Compound[i], " ", round(standards$lower[i], 5), " ", round(standards$upper[i], 5)))
  dev.off()
  
  allPeaks <- rbind(allPeaks, peaks3)
}
write.csv(allPeaks, "EICPeaks.csv", row.names = FALSE)

#### Bring back in EIC Peaks ####
EICPeaks <- readxl::read_excel("EICPeaks.xlsx")
EICPeaks %>% group_by(compound) %>% summarize(cv = sd(into) / mean(into))

########### All peaks ############
df0Peaks <- findChromPeaks(df0, CentWaveParam(ppm = 10, snthresh = 10, peakwidth = c(15, 40), prefilter = c(1, 5000)))
df0Spectra <- chromPeakSpectra(df0Peaks)
df0Peaks
chromPeaks(df0Peaks, mz = 216.1010, ppm = 5)
df0PeaksG <- groupChromPeaks(df0Peaks, param = PeakDensityParam(sampleGroups = c(1,1,1),
                                                                bw = 30, minFraction = .4, binSize = .1))
featureValues(df0PeaksG)
df0GSpectra <- chromPeakSpectra(df0PeaksG)

Atrazine1 <- df0Spectra[mcols(df0Spectra)$peak_id == "CP009079"]
plot(Atrazine1[[1]], Atrazine1[[1]])
AtrazineSepc <- data.frame(mz = Atrazine1[[1]]@mz, intensity = Atrazine1[[1]]@intensity)
ggplot(AtrazineSepc, aes(x = mz, xend = mz, y = 0, yend = intensity)) + geom_segment() + theme_bw()

## the map
ms1 <- which(msLevel(df1$`1`) == 1)
rtsel1 <- rtime(df1$`1`)[ms1] > 425 & rtime(df1$`1`)[ms1] < 480
M <- MSmap(object = df1$`1`, scans = ms1[rtsel1] , lowMz = 215, highMz =  217, resMz = .005, hd, zeroIsNA = TRUE)
plot(M, aspect = 1, allTicks = FALSE)
plot3D(M, rgl = TRUE)
plot3D(M, rgl = FALSE)

ms1[rtsel1][1]
ms1[rtsel1][length(ms1[rtsel1])]

M2 <- MSmap(object = df1$`1`, scans = 1507:1510 , lowMz = 0 , highMz =  400, resMz = .05, hd, zeroIsNA = TRUE)
plot3D(M2, rgl = FALSE)

ms2 <- which(msLevel(df1$`1`) == 2)
rtsel2 <- rtime(df1$`1`)[ms2] > 425 & rtime(df1$`1`)[ms2] < 480

precur2 <- precursorMz(df1$`1`[ms2[rtsel2]])
precur2 <- precur2 > 216.0999 & precur2 < 216.1021
precur2[precur2]

plot(df1$`1`[[1626]]) + theme_bw()

msp <- readLines('~/MoNA-export-MassBank.msp')
grepl("Name: ", msp)
