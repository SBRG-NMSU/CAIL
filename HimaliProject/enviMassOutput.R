########### Prereqs ###########
outputDir <- "C:/Users/ptrainor/Documents/GitHub/cail/HimaliProject/"
EMProjectDir <- "C:/Users/ptrainor/Dropbox (NMSU Advanced)/Patrick/Himali/enviMassProc/results/"
# Must run in R 3.3.3
library(enviMass)
setwd(outputDir)

########### Import and Process ###########
# Import profile data from enviMass results:
load_profileList(paste0(EMProjectDir, "profileList_pos"))
profileList_pos

profs <- enviMass::profiles_to_matrix(profileList_pos)
peaks <- profileList_pos$peaks
profInfo <- profileList_pos$index_prof

# Profile links:
load(paste0(EMProjectDir, "links_profiles_pos"))
getLinks <- function(x){
  return(paste(x$group$direct, collapse = ";"))
}

linksList <- list()
for(i in 1:length(links_profiles_pos)){
  if(length(links_profiles_pos[[i]]) > 1 || !is.na(links_profiles_pos[[i]])){
    linksList[[i]] <- data.frame(profID = names(links_profiles_pos)[i], links = getLinks(links_profiles_pos[[i]]))
  }
}
linksDF <- do.call("rbind", linksList)

# Internal standards:
load(paste0(EMProjectDir, "screening/res_IS_pos_screen"))
load(paste0(EMProjectDir, "screening/results_screen_IS_pos"))

getISTD <- function(x){
  x <- x[[1]]
  lNames <- names(x)
  tDF1 <- matrix(nrow = length(x$`m/z`), ncol = length(lNames))
  tDF1 <- as.data.frame(tDF1)
  names(tDF1) <- lNames
  tDF1$Peaks <- x$Peaks[,2]
  tDF1$score_1 <- x$score_1
  tDF1$score_2 <- x$score_2
  tDF1$`ppm deviation` <- x$`ppm deviation`
  tDF1$`RT deviation from mean` <- x$`RT deviation from mean`
  tDF1$`rescale factor` <- x$`rescale factor`
  tDF1$`m/z` <- x$`m/z`
  tDF1$Intensity <- x$Intensity
  tDF1$RT <- x$RT
  tDF1$file_ID <- x$file_ID
  return(tDF1)
}

iSTDLists <- list()
for(i in 1:length(res_IS_pos_screen)){
  temp1 <- do.call("rbind", lapply(res_IS_pos_screen[[i]], getISTD))
  iSTDLists[[i]] <- temp1
  if(!is.null(temp1)){
    iSTDLists[[i]]$iSTD_ID <- names(res_IS_pos_screen)[i]
  }
}
iSTDDF <- do.call("rbind", iSTDLists)

save(peaks, profs, profInfo, iSTDDF, linksDF, file = paste0(outputDir, "RData/enviMassOutput_20200826.RData"))
