############ Prereqs ############
options(stringsAsFactors = FALSE, scipen = 900)
oldPar <- par()
os <- Sys.info()["sysname"]
baseDir <- ifelse(os == "Windows", "C:/Users/ptrainor/gdrive/CAIL/WaterProject/AE_MS2_CmpID/", 
                  "~/gdrive/CAIL/WaterProject/AE_MS2_CmpID/")
setwd(baseDir)

library(tidyverse)

############ Import textfiles ############
# Annotation data from msp files:
load("mspInfo.RData")

# In Silico:
inSilico <- read.delim("inSilico_20200530.txt", header = TRUE, sep = "\t")
inSilico <- inSilico[, !names(inSilico) %in% c("File.path", "Title", "MS1.count", "X")]
colNames1 <- gsub("\\.1", "", colnames(inSilico)[grepl("\\.1", colnames(inSilico))])
colNames1b <- c("File.name", "MSMS.count", "PRECURSORMZ", "PRECURSORTYPE")

# Formula enumeration:
formEn <- read.delim("formulaEnumeration_20200530.txt", header = TRUE, sep = "\t")
formEn <- formEn[, !names(formEn) %in% c("File.path", "Title", "MS1.count", "X")]
colNames1c <- gsub("\\.1", "", colnames(formEn)[grepl("\\.1", colnames(formEn))])
colNames1d <- c("File.name", "MSMS.count", "PRECURSORMZ", "PRECURSORTYPE")

# Spectral matches:
spectra <- read.delim("spectralSearch_20200530.txt", header = TRUE, sep = "\t")
spectra <- spectra[, !names(spectra) %in% c("File.path", "Title", "MS1.count", "X")]
colNames1e <- gsub("\\.1", "", colnames(spectra)[grepl("\\.1", colnames(spectra))])
colNames1f <- c("File.name", "MSMS.count", "PRECURSORMZ", "PRECURSORTYPE")

############ Process in Silico ############
# In Silico wide to long:
sList1 <- c("", 1:4)
inSilicoTemp2 <- list()
for(i in 1:length(sList1)){
  if(i > 1){
    colNames2 <- paste0(colNames1, ".", sList1[i])
  }else{
    colNames2 <- colNames1
  }
  
  inSilicoTemp <- inSilico[, (names(inSilico) %in% colNames1b | names(inSilico) %in% colNames2)]
  inSilicoTemp$rank <- i
  names(inSilicoTemp) <- gsub("\\.\\d", "", names(inSilicoTemp))
  inSilicoTemp2[[i]] <- inSilicoTemp
}
inSilico <- do.call("rbind", inSilicoTemp2)
names(inSilico)[match(c("Structure.rank", "Total.score"), names(inSilico))] <- c("Structure", "StructureScore")
# Remove blank:
inSilico <- inSilico[!inSilico$Structure == "",]
# Cleanup duplicates:
inSilico <- inSilico %>% group_by(File.name, Structure, StructureScore, Formula, Ontology, SMILES) %>% 
  mutate(n1 = row_number(), n2 = n()) %>% ungroup()
inSilico$delete <- ifelse(inSilico$n2 > 1 & inSilico$n1 == 1, TRUE, FALSE)
inSilico <- inSilico %>% filter(delete == FALSE) %>% select(-n1, -n2, -delete)
rm(sList1, inSilicoTemp, inSilicoTemp2, colNames1, colNames1b, colNames2, i)

# Join with RT and precursor m/z:
df1$file <- gsub("\\.msp", "", df1$file)
inSilico <- df1 %>% left_join(inSilico, by = c("file" = "File.name"))

############ Process Formula enumeration ############
# Formula wide to long:
sList1 <- c("", 1:4)
formEnTemp2 <- list()
for(i in 1:length(sList1)){
  if(i > 1){
    colNames2 <- paste0(colNames1c, ".", sList1[i])
  }else{
    colNames2 <- colNames1c
  }
  
  formEnTemp <- formEn[, (names(formEn) %in% colNames1d | names(formEn) %in% colNames2)]
  formEnTemp$rank <- i
  names(formEnTemp) <- gsub("\\.\\d", "", names(formEnTemp))
  formEnTemp2[[i]] <- formEnTemp
}
formEn <- do.call("rbind", formEnTemp2)
names(formEn)[match(c("Formula.rank", "Formula.score"), names(formEn))] <- c("Formula", "FormulaScore")
formEn <- formEn[!formEn$Formula == "",]
rm(formEnTemp, formEnTemp2, colNames1d, colNames1c, colNames2, i, sList1)

# Join formula to inSilico
formEn <- formEn %>% select(file = File.name, formula = Formula, Theoretical.mass, massError = Mass.error, 
                         formulaScore = FormulaScore) %>% 
  left_join(inSilico %>% select(file, precursorMZ = precursormz, RT = rt, MSMSCount = MSMS.count, 
                           adduct = PRECURSORTYPE, structureScore = StructureScore, strucureRank = rank,
                           databases = Databases, formula = Formula, structure = Structure,
                           ontology = Ontology, InChIKey, SMILES),
            by = c("file", "formula"))

# Process for output:
formEn$RT <- as.numeric(gsub("RETENTIONTIME: ", "", formEn$RT))
formEn$precursorMZ <- as.numeric(gsub("PRECURSORMZ: ", "", formEn$precursorMZ))
formEn <- formEn %>% select(file, precursorMZ, RT, adduct, formula, theoreticalMass = Theoretical.mass, massError, 
                            formulaScore, structure, ontology, structureScore, InChIKey, SMILES, databases)
# Combined score:
formEn$combinedScore <- formEn$formulaScore + formEn$structureScore

# File name to feature number:
formEn$file <- as.integer(gsub("_Unknown", "", formEn$file))
names(formEn)[names(formEn) == "file"] <- "feature"
writexl::write_xlsx(formEn, "formEn_20200601.xlsx")

############ Process spectral matching ############
# Formula wide to long:
sList1 <- c("", 1:4)
spectraTemp2 <- list()
for(i in 1:length(sList1)){
  if(i > 1){
    colNames2 <- paste0(colNames1e, ".", sList1[i])
  }else{
    colNames2 <- colNames1e
  }
  
  spectraTemp <- spectra[, (names(spectra) %in% colNames1f | names(spectra) %in% colNames2)]
  spectraTemp$rank <- i
  names(spectraTemp) <- gsub("\\.\\d", "", names(spectraTemp))
  spectraTemp2[[i]] <- spectraTemp
}
spectra <- do.call("rbind", spectraTemp2)
names(spectra)[match(c("Structure.rank", "Total.score"), names(spectra))] <- c("Structure", "StructureScore")
spectra <- spectra[!spectra$Formula == "",]
rm(spectraTemp, spectraTemp2, colNames1e, colNames1f, colNames2, i, sList1)

# Process for export:
# Join with RT and precursor m/z:
spectra <- df1 %>% right_join(spectra, by = c("file" = "File.name"))
spectra <- spectra %>% select(file, precursorMZ = PRECURSORMZ, RT = rt, adduct = PRECURSORTYPE, Structure, StructureScore)
spectra$RT <- as.numeric(gsub("RETENTIONTIME: ", "", spectra$RT))
spectra <- spectra %>% filter(StructureScore > 0.65)

# File name to feature number:
spectra$file <- as.integer(gsub("_Unknown", "", spectra$file))
names(spectra)[names(spectra) == "file"] <- "feature"
writexl::write_xlsx(spectra, "spectra_20200601.xlsx")

############ Analysis of what is present ############
formEnTop <- formEn %>% arrange(feature, desc(combinedScore))
formEnTop <- formEnTop %>% group_by(feature) %>% filter(row_number() == 1)
ontologyDF <- as.data.frame(table(formEnTop$ontology))
names(ontologyDF) <- c("Ontology", "Count")
writexl::write_xlsx(ontologyDF, "ontologyDF.xlsx")
