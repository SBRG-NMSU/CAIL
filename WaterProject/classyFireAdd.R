############ Prereqs ############
# Begin always run
options(stringsAsFactors = FALSE, scipen = 600, max.print = 100000)
oldPar <- par()

library(tidyverse)

os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/", "~/gdrive/")
setwd(paste0(baseDir, "CAIL/WaterProject"))
# End always run

############ Get classyfire classification for top hits ############
# Import MetFrag results:
load(file = "RData/metFragCand_20200804.RData")

hitsList <- list()
for(i in 1:nrow(topHit)){
  if(topHit$InChI[i] != ""){
    hit <- classyfireR::get_classification(topHit$InChIKey[i])
    if(!is.null(hit)){
      hitClassDF <- hit@classification
      hitClassDF$Level <- factor(hitClassDF$Level, levels = hitClassDF$Level)
      hitClassDF <- hitClassDF %>% select(-CHEMONT) %>% spread(key = Level, value = Classification)
      hitClassDF$prof <- topHit$prof[i]
      hitsList[[i]] <- hitClassDF
      cat("Successful: ", i, "\n")
    }
    else{
      cat("Not successful: ", i, "\n")
    }
  }
}
hitList <- do.call("bind_rows", hitsList)
hitList <- hitList %>% select(prof, kingdom, superclass, class, subclass, level5 = `level 5`, 
                              level6 = `level 6`, level7 = `level 7`, level8 = `level 8`, level9 = `level 9`)
topHit <- topHit %>% left_join(hitList)

save(allHits, topHit, file = "RData/metFragCand_20200927.RData")
