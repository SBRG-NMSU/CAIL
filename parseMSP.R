options(stringsAsFactors = FALSE, scipen = 900)

library(tidyverse)
library(xcms)

msp <- readLines('~/MoNA-export-MassBank.msp')
bounds <- which(grepl("Name: ", msp))
mspList <- list()
bounds <- c(bounds, length(msp))
for(i in 1:(length(bounds) - 1)){
  mspList[[i]] <- msp[bounds[i]:(bounds[i + 1] - 3)]
}

mspList[[1]]

# Find where spectra stars:
end1 <- sapply(mspList, function(x) which(grepl("Num Peaks:", x)))
end2 <- sapply(mspList, function(x) length(x))

# Get just spectra:
mspListSpectra <- list()
for(i in 1:length(mspList)){
  mspListSpectra[[i]] <- mspList[[i]][(end1[i] + 1):end2[i]]
}

# Split spectra into m/z and Int
mspSpecDFFun <- function(x){
  specDF <- as.data.frame(str_split(mspListSpectra[[1]], " ", simplify = TRUE))
  specDF$V1 <- as.numeric(specDF$V1)
  specDF$V2 <- as.numeric(specDF$V2)
  names(specDF) <- c("mz", "int")
  return(specDF)
}
mspListSpectra2 <- lapply(mspListSpectra, mspSpecDFFun)
mspListSpectra <- mspListSpectra2
rm(mspListSpectra2)

names(mspListSpectra) <- paste0("spec", 1:length(mspListSpectra))
names(mspList) <- paste0("spec", 1:length(mspListSpectra))

annoDF <- data.frame(specID = paste0("spec", 1:length(mspListSpectra)))
annoDF$name <- annoDF$InChIKey <- annoDF$adduct <- ""
for(i in 1:nrow(annoDF)){
  annoDF$name[i] <- gsub("Name: ", "", mspList[[i]][grepl("Name: ", mspList[[i]])])
  annoDF$InChIKey[i] <- gsub("InChIKey: ", "", mspList[[i]][grepl("InChIKey: ", mspList[[i]])])
}
