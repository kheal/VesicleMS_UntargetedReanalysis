library(tidyverse)


#Set your datafiles
#Read in your output files from Skyline
Datfile1 <- "Intermediates/Targeted_QE_fromLauraC/HILIC-Pos_ProVesicles_QE.csv"
Datfile2 <- NA
Datfile3 <- NA
Datfile4 <- NA

#Read in a csv that has which samples go to which blanks
BlankMatcherFile <- "MetaDat/BlankMatcher.csv"

#Set your flag to find standards to base RT off of
StdFlag <- c("170124_Std_4uM_StdsInH2O_1", "170124_Std_4uM_StdsInMatrix_1", "170124_Std_4uM_StdsInMatrix_3", "170124_Std_4uM_StdsInH2O_2")

#Set your parameters
SNmin = 5
ppmflex = 6
Areamin = 40000  #40000 is good target for HILICNeg
RTflex = 1
BlankRatiomax = 10

#What's your output file name?
fileout <- "Intermediates/QCd_Tar_HILICPos.csv"
 
#Compounds to dump from this fraction
compounds.to.dump <- c("Adenosyl Homocysteine", "Adenosyl Methionine", "Glutathione")

#Says which are blanks and which are samples
BlankMatcher <- read_csv(BlankMatcherFile)

#Combine all the Datfiles
Datfile_Comb <- read.csv(Datfile1) %>% 
  unique()

#First do easy flags - SN, ppm, Areamin
Datorig <- Datfile_Comb %>%
  select(-Protein.Name, -Protein) %>%
  mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
  mutate(Area = as.numeric(as.character(Area))) %>%
  mutate(Background = as.numeric(as.character(Background))) %>%
  mutate(Mass.Error.PPM = as.numeric(as.character(Mass.Error.PPM))) 

Dat1 <- Datorig %>%
  filter(Replicate.Name %in% BlankMatcher$Replicate.Name) %>%
  mutate(SNFlag = ifelse((Area/Background < SNmin), "SNFlag", NA)) %>%
  mutate(ppmFlag = ifelse((abs(Mass.Error.PPM) > ppmflex), "ppmFlag", NA)) %>%
  mutate(areaminFlag = ifelse((Area < Areamin), "areaminFlag", NA))

#Next do RT flag compared to standard or a template
RT_of_Stds <- Datorig %>%
  filter(Replicate.Name %in% StdFlag) %>%
  select(Precursor.Ion.Name, Retention.Time) %>%
  group_by(Precursor.Ion.Name) %>%
  summarise(RT_ref = mean((Retention.Time), na.rm = TRUE))

Dat2 <- Dat1 %>%
  left_join(RT_of_Stds, by = "Precursor.Ion.Name") %>%
  mutate(Retention.Time = as.numeric(as.character(Retention.Time))) %>%
  mutate(RTFlag = ifelse((abs(Retention.Time - RT_ref) > RTflex), "RTFlag", NA))

#Next do Blank flag compared to blank
Area_of_Blanks <- Datorig %>%
  filter(Replicate.Name %in% BlankMatcher$Blank.Name) %>%
  rename(Blank.Name = Replicate.Name,
         Blank.Area = Area) %>%
  select(Blank.Name, Precursor.Ion.Name, Blank.Area) %>%
  left_join(BlankMatcher, by = "Blank.Name") %>% select(-Blank.Name) %>%
  arrange(desc(Blank.Area)) %>%
  group_by(Precursor.Ion.Name, Replicate.Name) %>% filter(row_number() == 1)

Dat3 <- Dat2 %>%
  left_join(Area_of_Blanks, by = c("Replicate.Name", "Precursor.Ion.Name")) %>%
  mutate(BlankFlag = ifelse(Area/Blank.Area < BlankRatiomax, "BlankFlag", NA))

#Finally, combine all the flags and throw out any peak with a flag
Dat4 <- Dat3 %>%
  mutate(Flags = paste(SNFlag, ppmFlag, areaminFlag, RTFlag, BlankFlag, sep = ", ")) %>%
  mutate(Flags = as.character(Flags %>% str_remove_all("NA, ") %>%  str_remove_all("NA"))) %>%
  mutate(QC_area = ifelse(str_detect(Flags, "Flag"), NA, Area))%>%
  filter(!Precursor.Ion.Name %in% compounds.to.dump)

#To inspect
Dat5 <- Dat4 %>%
  select(Precursor.Ion.Name, Replicate.Name, Flags, QC_area)

#To write the file
comment.text <- paste("Hello! Welcome to the world of QE Quality Control! ",
                      "Minimum area for a real peak: ", Areamin, ". ",
                      "RT flexibility: ", RTflex, ". ",
                      "Blank can be this fraction of a sample: ", BlankRatiomax, ". ",
                      "S/N ratio: " , SNmin, ". ",
                      "Parts per million flexibility: ", ppmflex, ". ",
                      "Processed on: ", Sys.time(), ". ",
                      sep = "")
new.filename <- fileout
con <- file(new.filename, open="wt")
writeLines(paste(comment.text), con)
write.csv(Dat4, con)
close(con)
