options(stringsAsFactors = FALSE, scipen = 900)

setwd("C:/Users/ptrainor/gdrive/CAIL/WaterProject")

msp <- readLines('AE_MS2_ID/Algal_eff_All-27samples_50mlSPE_DFX75_Chathurica_MS2_1.msp')
bounds <- which(grepl("(Name|NAME)", msp))

mspList <- list()
bounds <- c(bounds, length(msp))
for(i in 1:(length(bounds) - 1)){
  mspList[[i]] <- msp[bounds[i]:(bounds[i + 1] - 1)]
}

# Throw out ones w/o MS2 data:
whichHaveMS2 <- which(sapply(mspList, function(x) sum(grepl("Num Peaks: 0", x))) == 0)
mspList2 <- mspList[whichHaveMS2]

mspList2 <- do.call("c", mspList2)

write.table(mspList2, file = 'AE_MS2_ID/Algal_eff_All-27samples_50mlSPE_DFX75_Chathurica_MS2_1_v2.msp', sep = "\n", quote = FALSE,
            row.names = FALSE, col.names = FALSE, eol = "\n")
