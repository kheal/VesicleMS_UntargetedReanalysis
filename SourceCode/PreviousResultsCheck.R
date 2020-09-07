library(ggplot2)
library(tidyverse)
library(fuzzyjoin)
options(dplyr.summarise.inform=F)
options(readr.num_columns = 0)

#Define all your inputs here------
Laura.results.file <- "Intermediates/From_LauraH/AllMFs_Ranks_Final_withNotes.csv"
untargeted.results <- "Intermediates/FullMF_biovolNormalized.csv"


#Check out that dat
Laura.results.full <- read_csv(Laura.results.file) %>%
  filter(To_Check == "Y") %>%
  mutate(fraction.9312 = fraction.9312 %>% str_replace("CyanoDCM", "CyanoOrg"))
my.results.full <- read_csv(untargeted.results)

#See if there are any matches - first do it for just the CyanoDCM, dump RT for now
fractions <- c("HILICNeg", "HILICPos", "CyanoOrg", "CyanoAq")
overlap.test <- list()

for (i in 1:4){
Laura.results.sub <- Laura.results.full %>%
  select(MassFeature.9312, mz.9312, RT.9312, fraction.9312, Compound.Name.9312) %>%
  rename(Fraction = fraction.9312, 
         mz = mz.9312, 
         Retention_time = RT.9312) %>% 
  filter(Fraction == fractions[i]) %>%
  select(-Fraction)

my.results.sub <- my.results.full %>%
  select(MF, Fraction, mz, Retention_time) %>%
  filter(Fraction == fractions[i]) %>%
  select(-Fraction)

overlap.test[[i]] <- my.results.sub %>%
  difference_inner_join(Laura.results.sub, by ="mz", max_dist = 0.005) %>%
  mutate(Fraction = fractions[i]) %>%
  mutate(RTdiff = Retention_time.x-Retention_time.y) %>%
  mutate(lowRT = ifelse(Fraction %in% fractions[1:2] & abs(RTdiff) < 0.5, TRUE,
                        ifelse(Fraction %in% fractions[3:4] & abs(RTdiff) < 0.008, TRUE, FALSE)))%>%
  filter(lowRT == TRUE) %>%
  mutate(Fraction == fractions[i])
}
overlap.test.df <- do.call(bind_rows, overlap.test)
overlap.MFs <- overlap.test.df %>%
  select(Fraction, MF, RTdiff, Retention_time.x, mz.x, Compound.Name.9312)%>%
  unique()

#Check if they are present in the vesicles
vesicle.check <- my.results.full %>%
  filter(MF %in% overlap.MFs$MF)
