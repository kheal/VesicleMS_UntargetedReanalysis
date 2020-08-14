library(ggplot2)
library(tidyverse)
options(dplyr.summarise.inform=F)
options(readr.num_columns = 0)

#Define all your inputs here------
samp.key.file <- "MetaDat/SampleKey_Vesicles_QE.csv"
HILIC.longdat.file <- "Intermediates/BMISReports/BMISd_Areas_long_HILIC.csv"
CyanoAq.longdat.file  <- "Intermediates/BMISReports/BMISd_Areas_long_Cyano.csv"
CyanoOrg.longdat.file <- "Intermediates/BMISReports/BMISd_Areas_long_CyanoOrg.csv"
smp.to.blank.ratio <- 10
rsd.of.poo.max <- 0.3
min.peak.area <- 5000
min.SN <- 5
HILIC.SN.file.1 <- "HILICPos/SN_0_20208111327.txt"
HILIC.SN.file.2 <- "HILICNeg/SN_0_20208111351.txt"
CyanoAq.SN.file <- "RP/SN_0_20208111421.txt"
CyanoOrg.SN.file <- "RP_org/SN_2_20208131211.txt"
cell.count.file <- "MetaDat/cell_counts.csv"

#Load and clean up the SNs-----
SN.dat.HILIC.1 <- read_delim(HILIC.SN.file.1,
                           "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "HILICPos") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`) %>%
  select(MF, Retention_time, mz, `170124_Blk_ProcessBlk_1`:`170124_Smp_Vesicle9313_3`) %>%
  pivot_longer(`170124_Blk_ProcessBlk_1`:`170124_Smp_Vesicle9313_3`, names_to = "SampID", values_to = "SN") %>% 
  mutate(SampID = SampID %>% 
           str_replace("_FS_","_FS") %>%
           str_replace("Jan25_","Jan25") %>%
           str_replace("Jan24_","Jan24") %>%
           str_replace("_1a", "_vesicleblank_1") %>%
           str_replace("_2a", "_vesicleblank_2") %>%
           str_replace("_3a", "_vesicleblank_3") %>%
           str_replace("_4a", "_vesicleblank_4")) 
SN.dat.HILIC.2 <- read_delim(HILIC.SN.file.2,
                             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "HILICNeg") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`) %>%
  select(MF, Retention_time, mz,`170124_Blk_ProcessBlk_1`:`170124_Smp_Vesicle9313_3`) %>%
  pivot_longer(`170124_Blk_ProcessBlk_1`:`170124_Smp_Vesicle9313_3`, names_to = "SampID", values_to = "SN") %>% 
  mutate(SampID = SampID %>% 
           str_replace("_FS_","_FS") %>%
           str_replace("Jan25_","Jan25") %>%
           str_replace("Jan24_","Jan24") %>%
           str_replace("_1a", "_vesicleblank_1") %>%
           str_replace("_2a", "_vesicleblank_2") %>%
           str_replace("_3a", "_vesicleblank_3") %>%
           str_replace("_4a", "_vesicleblank_4")) 
SN.dat.CyanoAq <- read_delim(CyanoAq.SN.file,
                             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "CyanoAq") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`) %>%
  select(MF, Retention_time, mz, `170128_Blk_1a`:`170128_Smp_Vesicle9313_3`) %>%
  pivot_longer(`170128_Blk_1a`:`170128_Smp_Vesicle9313_3`, names_to = "SampID", values_to = "SN") %>% 
  mutate(SampID = SampID %>% 
           str_replace("_FS_","_FS") %>%
           str_replace("Jan25_","Jan25") %>%
           str_replace("Jan24_","Jan24") %>%
           str_replace("_1a", "_vesicleblank_1") %>%
           str_replace("_2a", "_vesicleblank_2") %>%
           str_replace("_3a", "_vesicleblank_3") %>%
           str_replace("_4a", "_vesicleblank_4"))
SN.dat.CyanoOrg  <- read_delim(CyanoOrg.SN.file,
                             "\t", escape_double = FALSE, trim_ws = TRUE,  skip = 4) %>%
  mutate(Column = "CyanoOrg") %>% 
  mutate(MF = paste0(`Alignment ID`, "_", Column)) %>%
  rename(Retention_time = `Average Rt(min)`,
         mz = `Average Mz`) %>%
  select(MF, Retention_time, mz, `170128_Blk_1a_DCM`:`170128_Smp_Vesicle9313_3_DCM`) %>%
  pivot_longer(`170128_Blk_1a_DCM`:`170128_Smp_Vesicle9313_3_DCM`, names_to = "SampID", values_to = "SN") %>% 
  mutate(SampID = SampID %>% 
           str_replace("_DCM","") %>%
           str_replace("_FS_","_FS") %>%
           str_replace("Jan25_","Jan25") %>%
           str_replace("Jan24_","Jan24") %>%
           str_replace("_1a", "_vesicleblank_1") %>%
           str_replace("_2a", "_vesicleblank_2") %>%
           str_replace("_3a", "_vesicleblank_3") %>%
           str_replace("_4a", "_vesicleblank_4"))

sn.dat.combo <-SN.dat.HILIC.1 %>% bind_rows(SN.dat.HILIC.2) %>%
  bind_rows(SN.dat.CyanoAq) %>% bind_rows(SN.dat.CyanoOrg)

MF.info.combo <-sn.dat.combo %>%
  select(MF, Retention_time, mz) %>% unique()


#Load your BMISd longdats match by samples, attach the SN data------
dat <- read_csv(HILIC.longdat.file) %>%
  bind_rows(read_csv(CyanoAq.longdat.file))%>%
  bind_rows(read_csv(CyanoOrg.longdat.file)) %>%
  left_join(sn.dat.combo, by = c("MF", "SampID"))

#For each MF, get the following information: sample type blank area; CV of pooled (per sample type)-----
dat.RSD.poo <- dat %>%
  filter(type == "Poo") %>%
  group_by(MF, Samp.Type) %>%
  summarise(RSD_poo = sd(Adjusted_Area)/mean(Adjusted_Area)) 

dat.blank.ave <- dat %>%
  filter(type == "Blk") %>%
  group_by(MF, Samp.Type) %>%
  summarise(ave_blank = mean(Adjusted_Area)) %>%
  mutate(ave_blank = ifelse(ave_blank == 0, 100, ave_blank))

dat.smp.ave <- dat %>%
  filter(type == "Smp") %>%
  group_by(MF, Samp.Type) %>%
  summarise(ave_smp = mean(Adjusted_Area)) 

dat.MF.good <- dat.RSD.poo %>%
  left_join(dat.blank.ave, by = c("MF", "Samp.Type")) %>%
  left_join(dat.smp.ave, by = c("MF", "Samp.Type")) %>%
  mutate(smptoblank = ave_smp/ave_blank) %>%
  mutate(GoodMFs = ifelse(smptoblank > smp.to.blank.ratio & RSD_poo < 0.3, 1, 0)) %>%
  mutate(GoodMFs = ifelse(is.na(GoodMFs), 0, GoodMFs))

#See if each individual peaks are good enough by minimum peak size and maximum signal to noise-----
dat.QC <- dat %>%
  mutate(GoodPeak = ifelse(Adjusted_Area > min.peak.area & SN > 5, 1, 0)) %>%
  left_join(dat.MF.good, by = c("MF", "Samp.Type"))


#Summarize number of "good mass features in each fraction and in each sample type"-----
dat.QC.summary <- dat.QC %>%
  mutate(Fraction = str_extract(MF, "_\\w+$") %>%
           str_replace("_", "")) %>%
  filter(type == "Smp") %>%
  select(MF, SampID:GoodPeak, GoodMFs, Fraction) 
dat.QC.summary.2 <- dat.QC.summary %>%
  group_by(MF, samp, Samp.Type, Fraction) %>%
  summarise(GoodPeakTotal = sum(GoodPeak),
            GoodMFsTotal = sum(GoodMFs)) %>%
  mutate(Keep = ifelse(GoodPeakTotal > 1.5 & GoodMFsTotal == 3, 1, 0)) 
test <- dat.QC.summary.2 %>% filter(is.na(Keep))

dat.QC.summary.3 <- dat.QC.summary.2 %>%
  group_by(samp, Samp.Type, Fraction) %>%
  summarise(MFs_observed = sum(Keep))
write_csv(dat.QC.summary.3, "Intermediates/MF_observed_summary.csv")


#Summarize number of "Shared between cells or vesicles of both strains"-----
dat.QC.summary.5 <- dat.QC.summary.2 %>%
  mutate(Strain = str_extract(samp, "931.")) %>%
  group_by(MF, Samp.Type, Fraction) %>%
  summarise(Keep_total = sum(Keep))%>%
  mutate(Shared = ifelse(Keep_total > 1.5, 1, 0)) %>%
  mutate(At_least_one = ifelse(Keep_total >0, 1, 0)) %>%
  group_by(Samp.Type, Fraction) %>%
  summarise(Shared_total_btwStrains = sum(Shared),
            Shared_percent_btwStrains = sum(Shared)/sum(At_least_one)*100)
write_csv(dat.QC.summary.5, "Intermediates/MF_observed_summary_betweenstrains.csv")


#Summarize number of "Shared between cells and vesicles of both strains"-----
dat.QC.summary.6 <- dat.QC.summary.2 %>%
  mutate(Strain = str_extract(samp, "931.")) %>%
  group_by(MF, Strain, Fraction) %>%
  summarise(Keep_total = sum(Keep)) %>%
  mutate(Shared = ifelse(Keep_total > 1.5, 1, 0)) %>%
  mutate(At_least_one = ifelse(Keep_total >0, 1, 0))%>%
  group_by(Strain, Fraction) %>%
  summarise(Shared_total_btwVandC = sum(Shared),
            Shared_percent_btwVandC = sum(Shared)/sum(At_least_one)*100)
write_csv(dat.QC.summary.6, "Intermediates/MF_observed_summary_betweensVandC.csv")

#Make supplemental table
dat.QC.2 <- dat.QC %>% left_join(read_csv(cell.count.file), by = c("samp", "replicate")) %>%
  mutate(Strain = str_extract(samp, "931.")) %>%
  filter(type == "Smp") %>%
  mutate(QCd_biovolume_area = ifelse(GoodPeak == 1 & GoodMFs == 1, Adjusted_Area/biovolume_extracted, NA))
dat.QC.2.wide <- dat.QC.2 %>%
  mutate(ColumnID = paste(Samp.Type, Strain, replicate, sep = "_")) %>%
  mutate(ColumnID = ColumnID %>% str_replace("Pellet", "Cell")) %>%
  select(MF, ColumnID, QCd_biovolume_area) %>%
  filter(!is.na(QCd_biovolume_area)) %>%
  pivot_wider(names_from = ColumnID, values_from = QCd_biovolume_area) %>%
  left_join(MF.info.combo) %>%
  mutate(Fraction = str_extract(MF, "_\\w+$") %>%
                                       str_replace("_", "")) %>%
  select(MF, Fraction, mz, Retention_time, everything()) %>%
  arrange(Fraction, mz)
write_csv(dat.QC.summary.6, "Intermediates/FullMF_biovolNormalized.csv")

