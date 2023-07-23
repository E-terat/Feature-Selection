#In this document we will construct the count matrix of all the dataset download from GEO.
#including HC (Hepatic), Cirr (cirrhosis), MGUS, MM (myeloma multiple) and HD (human donor) 
#samples

library(tidyverse)
source(here("utils/build_counts_matrix.R"))

#all metadata
metadata <- read.delim2("./data/metadata/GSE182824_series_matrix_v2.txt") %>% t()
metadata <- cbind(rownames(metadata), metadata)
colnames(metadata) <- metadata[1,]
metadata <- metadata[-1,]
rownames(metadata) <- NULL
metadata <- as.data.frame(metadata)
metadata$outcome<- as.factor(metadata$outcome) %>% 
  factor(labels = c("HC", "Cirr", "MGUS", "MM", "HD"))

##all counts_dedup with the results of the rna-seq pipeline explained in 'Material y Métodos' 
##in the main text
procesado_1 <- here("data/procesado_1/counts_dedup/")
procesado_2 <- here("data/procesado_2/counts_dedup/")
procesado_3 <- here("data/procesado_3/counts_dedup/")
procesado_4 <- here("data/procesado_4/counts_dedup/")
procesado_5 <- here("data/procesado_5/counts_dedup/")

dedupCounts1 <- build_counts_matrix(procesado_1)
dedupCounts2 <- build_counts_matrix(procesado_2)
dedupCounts3 <- build_counts_matrix(procesado_3)
dedupCounts3_re <- build_counts_matrix(procesado_3, samplePattern = "RSQ_...") 
col <- c(colnames(dedupCounts3), colnames(dedupCounts3_re))
col <- unique(col[!is.na(col)])
colnames(dedupCounts3) <- col  
dedupCounts4 <- build_counts_matrix(procesado_4, samplePattern = "RSQ_...")
dedupCounts5 <- build_counts_matrix(procesado_5, samplePattern = "RSQ_...")

# Construir la matriz de conteos global con todas las muestras disponibles
dedupCounts <- dedupCounts1 %>%
  full_join(dedupCounts2, by="GeneID") %>%
  full_join(dedupCounts3, by="GeneID") %>%
  full_join(dedupCounts4, by="GeneID") %>%
  full_join(dedupCounts5, by="GeneID") 


dir.create(here("data/dedupCounts_data"))
write_csv(dedupCounts, file = here("data/dedupCounts_data/dedupCounts_all.csv"))




rm(list=c("dedupCounts1", "dedupCounts2", "dedupCounts3", "dedupCounts3_re","procesado_1", "procesado_2", "procesado_3"))
