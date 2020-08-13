library(ggplot2)
library(tidyverse)
options(dplyr.summarise.inform=F)
options(readr.num_columns = 0)

#Define all your inputs here
samp.key.file <- "MetaDat/SampleKey_Vesicles_QE.csv"
HILIC.longdat.file <- "Intermediates/BMISReports/BMISd_Areas_long_HILIC.csv"
CyanoAq.longdat.file  <- "Intermediates/BMISReports/BMISd_Areas_long_Cyano.csv"
CyanoDCM.longdat.file <- "Intermediates/BMISReports/BMISd_Areas_long_CyanoOrg.csv"


#Load your files match by samples
dat <- read_csv(HILIC.longdat.file) %>%
  bind_rows(read_csv(CyanoAq.longdat.file)) %>%
  bind_rows(read_csv(CyanoDCM.longdat.file))

#For each MF, get the following information: sample type blank area; CV of pooled (per sample type)
dat.RSD.poo <- dat %>%
  filter(type == "Poo") %>%
  group_by(MF, Samp.Type) %>%
  summarise(RSD = sd(Adjusted_Area)/mean(Adjusted_Area))
