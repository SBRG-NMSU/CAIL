############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/CAIL/YuklProt/PrelimData/", 
                  "~/gdrive/CAIL/YuklProt/PrelimData/")
setwd(baseDir)

library(tidyverse)

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5))

############ Import and process text files ############
peptides <- readxl::read_excel("txt/peptides.xlsx", guess_max = 100000)
evidence <- readxl::read_excel("txt/evidence.xlsx", guess_max = 100000)
protGroups <- readxl::read_excel("txt/proteinGroups.xlsx", guess_max = 100000)
msmsScan <- readxl::read_excel("txt/msmsScans.xlsx", guess_max = 100000)
msms <- readxl::read_excel("txt/msms.xlsx", guess_max = 100000)
pSTY <- readxl::read_excel("txt/Phospho (STY)Sites.xlsx", guess_max = 100000)
modPep <- readxl::read_excel("txt/modificationSpecificPeptides.xlsx", guess_max = 100000)

swissprot <- read.csv("uniprot-reviewed yes+taxonomy 318586.csv", sep="\t")

############ Possible phosphorylated protein ############
poss1 <- protGroups %>% filter(id == 364)
poss1Pep <- peptides %>% filter(id %in% str_split(poss1$`Peptide IDs`, ";", simplify = TRUE))
poss1ModPep <- modPep %>% filter(id %in% str_split(poss1$`Mod. peptide IDs`, ";", simplify = TRUE))

rm(poss1, poss1Pep, poss1ModPep)

############ Protein processing and normalization ############
length(unique(protGroups[str_split(protGroups$`Protein IDs`, ";", simplify = TRUE)[,2] == "",]$`Protein IDs`)) # 1,243
allProts <- unique(c(as.matrix(str_split(protGroups$`Protein IDs`, ";", simplify = TRUE))))
table(swissprot$Entry %in% allProts)
prop.table(table(swissprot$Entry %in% allProts))

# Add present vs absent:
protGroups$presWT <- ifelse(protGroups$`Intensity WT` > 0, TRUE, FALSE)
protGroups$presHNOX <- ifelse(protGroups$`Intensity HNOX` > 0, TRUE, FALSE)
xtabs(~presWT + presHNOX, data = protGroups)

# Normalization and FC for present in both:
protGroups2 <- protGroups %>% filter(presWT & presHNOX) %>% select(id, `Intensity WT`, `Intensity HNOX`)
protVsn <- vsn::justvsn(as.matrix(protGroups2[, 2:3]))
protGroups2$IntWT <- protVsn[,"Intensity WT"]
protGroups2$IntHNOX <- protVsn[,"Intensity HNOX"]
protGroups2$FC <- protGroups2$IntHNOX - protGroups2$IntWT
protGroups2$FCType <- ifelse(protGroups2$FC > 1, "Higher", ifelse(protGroups2$FC < -1, "Lower", "Ind"))
protFC <- xtabs(~ FCType, data = protGroups2)

# Highest FC:
protGroups3 <- protGroups2 %>% arrange(desc(abs(FC)))
protGroups3 <- protGroups3[1:5,]
protGroups4 <- protGroups3 %>% left_join(protGroups[protGroups$id %in% protGroups3$id,])
protGroups4 <- protGroups4 %>% arrange(FC)

# Get peptides for those:
pG3pep <- peptides[peptides$id %in% do.call("c", sapply(protGroups[protGroups$id %in% protGroups3$id,]$`Peptide IDs`, function(x) str_split(x, ";"))),]
pG3pepE <- evidence[evidence$`Peptide ID` %in% do.call("c", sapply(protGroups[protGroups$id %in% protGroups3$id,]$`Peptide IDs`, function(x) str_split(x, ";"))),]

# Protein groups vs proteins:
length(unique(protGroups[str_split(protGroups$`Protein IDs`, ";", simplify = TRUE)[,2] == "",]$`Protein IDs`))

# All protein DE:
protGroups2b <- protGroups %>% select(id, `Intensity WT`, `Intensity HNOX`, `Protein IDs`,
                                      `Fasta headers`)
protVsnb <- vsn::justvsn(as.matrix(protGroups2b[, 2:3]))
protGroups2b$intWT <- protVsnb[,"Intensity WT"]
protGroups2b$intHNOX <- protVsnb[,"Intensity HNOX"]
protGroups2b$FC <- protGroups2b$intHNOX - protGroups2b$intWT
writexl::write_xlsx(protGroups2b, path = "Contrasts_20200910.xlsx")

############ Peptide processing and normalization ############
# Present vs. absent:
peptides2 <- peptides %>% select(id, `Intensity WT`, `Intensity HNOX`) %>% as.data.frame()
peptides2$presWT <- ifelse(peptides2$`Intensity WT` > 0, TRUE, FALSE)
peptides2$presHNOX <- ifelse(peptides2$`Intensity HNOX` > 0, TRUE, FALSE)
peptidesPresAbs <- xtabs(~ presHNOX + presWT, data = peptides2)

# Normalization and FC for present in both:
peptides2 <- peptides2 %>% filter(presWT & presHNOX)
pepVsn <- vsn::justvsn(as.matrix(peptides2[, 2:3]))
peptides2$IntWT <- pepVsn[,"Intensity WT"]
peptides2$IntHNOX <- pepVsn[,"Intensity HNOX"]
peptides2$FC <- peptides2$IntHNOX - peptides2$IntWT
peptides2$FCType <- ifelse(peptides2$FC > 1, "Higher", ifelse(peptides2$FC < -1, "Lower", "Ind"))
peptidesFC <- xtabs(~ FCType, data = peptides2)
