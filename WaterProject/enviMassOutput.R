library(enviMass)
setwd("C:/Users/ptrainor/gdrive/CAIL/WaterProject")

# Import profile data from enviMass results:
load_profileList("C:/Users/ptrainor/gdrive/CAIL/WaterProject/WaterProjectDemo2/results/profileList_pos")
profileList_poss

profs <- enviMass::profiles_to_matrix(profileList_pos)
peaks <- profileList_pos$peaks
# peaks3753 <- peaks[peaks[,"profileIDs"] == 3753, ]
