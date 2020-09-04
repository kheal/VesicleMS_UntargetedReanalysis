library(ggplot2)
library(tidyverse)
options(dplyr.summarise.inform=F)
options(readr.num_columns = 0)

#Define all your inputs here------
samp.key.file <- "MetaDat/SampleKey_Vesicles_QE.csv"
cell.count.file <- "MetaDat/cell_counts.csv"
HILICPos.target.file <- "Intermediates/QCd_Tar_HILICPos.csv"
HILICNeg.target.file <- "Intermediates/QCd_Tar_HILICNeg.csv"
RPAq.target.file <- "Intermediates/QCd_Tar_RPAq.csv"
RPOrg.target.file <- "Intermediates/QCd_Tar_RPOrg.csv"

#Read in Samp Key
samp.key <- read_csv(samp.key.file)
cell.count <- read_csv(cell.count.file)

#Combine the datasets and add fraction information
combo.dat <- read_csv(HILICPos.target.file, skip = 1) %>% mutate(Fraction = "HILICPos") %>% select(-Height) %>%
  bind_rows(read_csv(HILICNeg.target.file, skip = 1) %>% mutate(Fraction = "HILICNeg")%>% select(-Height))%>%
  bind_rows(read_csv(RPAq.target.file, skip = 1) %>% mutate(Fraction = "RPAq")%>% select(-Height))%>%
  bind_rows(read_csv(RPOrg.target.file, skip = 1) %>% mutate(Fraction = "RPOrg")%>% select(-Height))


#Clean up combo dat
combo.dat.2 <- combo.dat %>%
  rename(sample = Replicate.Name,
         compound = Precursor.Ion.Name) %>%
  select(sample, compound, Fraction, QC_area) %>%
  left_join(samp.key %>% rename(sample = Sample.Name) %>% select(sample, Samp.Type), by = "sample")%>%
  mutate(cell.type = str_extract(sample, "931\\d"))

#Add cell count info
cell.count.sub <- cell.count %>% 
  select(-samp) %>% 
  mutate(cell.type = as.character(cell.type)) %>%
  unique()

combo.dat.3 <- combo.dat.2 %>%
  separate(sample, into= c("date", "Smp", "sampID", "replicate"), remove = FALSE) %>%
  select(-date, -Smp, -sampID) %>%
  left_join(cell.count.sub) %>%
  mutate(peakarea_perCell = QC_area/biovolume_extracted)
  
#Summarize data
summary.dat <- combo.dat.3 %>%
  mutate(observed = ifelse(peakarea_perCell > 0, 1, 0)) %>%
  group_by(compound, Samp.Type, Fraction, cell.type) %>%
  summarise(AvePeakArea = mean(peakarea_perCell, na.rm = TRUE),
           StdDevPeakArea = sd(peakarea_perCell, na.rm = TRUE),
           numObserved = sum(observed, na.rm = TRUE)) 

#get list of unobserved compounds
never.seen.compounds <- combo.dat.3 %>%
  mutate(observed = ifelse(peakarea_perCell > 0, 1, 0)) %>%
  group_by(compound, Fraction) %>%
  summarise(numObserved = sum(observed, na.rm = TRUE)) %>%
  filter(numObserved == 0)

#Make it like Supp Table 7
summary.dat.2 <- summary.dat %>%
  filter(!compound %in% never.seen.compounds$compound)%>%
  mutate(col.id = paste(cell.type, Samp.Type, sep = "_"))%>%
  mutate(summary_stat = paste0(signif(AvePeakArea, digits = 2), " (", signif(StdDevPeakArea, digits = 2), "), n=", numObserved)) %>%
  mutate(summary_stat = ifelse(AvePeakArea>0, summary_stat, NA)) %>%
  ungroup() %>%
  select(compound, Fraction, col.id, summary_stat) %>%
  pivot_wider(names_from = col.id, values_from = summary_stat)


write_csv(summary.dat.2, "Intermediates/Targeted_summary.csv")
  