library(tidyverse)
library(here)

#In this document we will construct the metadata file for the Myeloma Multiple dataset 

#as the main dataset contained Hepatocellular cancer samples, the code below includes
#the separation of both myeloma and hepatic datasets

## metadata 
metadata <- utils::read.csv(here("data/metadata/all_metadata.csv"))  %>%
  mutate(sample = str_remove(sample, "b"), 
         outcome = as.factor(outcome), 
         age = as.numeric(age), 
         gender = as.factor(gender),
         RNA_extraction = as.factor(RNA_extraction),
         Library_preparation = as.factor(Library_preparation), 
         set = as.factor(set), 
         race = as.factor(race),
         group = as.factor(group)) %>%
  filter(sample != "PP084")   #no information of this sample


metadata <- metadata %>% 
  filter(sample != "PP003",   #races != WHITE
         sample != "PP026",
         sample != "PP029",
         sample != "RSQ_047")

write.csv(metadata, "./data/metadata/all_metadata_prep.csv")


## construct the count matrix of the myeloma dataset
MM_samples <- metadata %>% filter(outcome == "HD" | outcome == "MM")
HC_samples <- metadata %>% filter(outcome == "HD" | outcome == "HC")

dedupCounts_MM <- utils::read.csv("./data/dedupCounts_data/dedupCounts_all.csv") %>% 
  select(GeneID, MM_samples$sample) 
write_csv(dedupCounts_MM, "./data/dedupCounts_data/dedupCounts_mm.csv")

dedupCounts_HC <- utils::read.csv("./data/dedupCounts_data/dedupCounts_all.csv") %>%
  select(GeneID, HC_samples$sample) 
write_csv(dedupCounts_HC, "./data/dedupCounts_data/dedupCounts_hc.csv")

