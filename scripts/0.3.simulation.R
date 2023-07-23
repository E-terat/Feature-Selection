
library(here)
library(tidyverse)
### genes difexpressed -------------------------------------
set.seed(12)
size <- 200
prob <- 0.99

n_genes <- 5
n_replicas <- 20

simulated_data_control <- matrix(data = rnbinom(n_genes * n_replicas, size = size, prob = prob), ncol = n_replicas)

size <- 200
prob <- 0.8

n_genes <- 5
n_replicas <- 20

simulated_data_condicion <- matrix(data = rnbinom(n_genes * n_replicas, size = size, prob = prob), ncol = n_replicas)

simulated_target_genes <- cbind(simulated_data_control, simulated_data_condicion)
simulated_target_genes[4,] <- simulated_target_genes[4,][40:1]
simulated_target_genes[5,] <- simulated_target_genes[5,][40:1]

# gene names 
gene_names <- paste0("Gene", 1:n_genes)

# sample names 
replica_names <- paste0("Sample", 1:(n_replicas*2))

rownames(simulated_target_genes) <- gene_names
colnames(simulated_target_genes) <- replica_names


### non diffexpressed genes----------------------------------
size <- 100
prob <- 0.95

n_genes <- 100
n_replicas <- 40

simulated_data2 <- matrix(data = rnbinom(n_genes * n_replicas, size = size, prob = prob), 
                          ncol = n_replicas)

# gene names
gene_names <- paste0("Gene", 6:(n_genes + 5))

# sample names
replica_names <- paste0("Sample", 1:n_replicas)

# Asignar nombres a las filas y columnas de la matriz
rownames(simulated_data2) <- (gene_names)
colnames(simulated_data2) <- replica_names

### all dataset -------------------------------------------

all_dataset <- rbind(simulated_target_genes, simulated_data2) %>% 
  data.frame() %>%
  rownames_to_column(var = "GeneID")


# save data
write_csv(all_dataset, file = here("data/simulated_data/sim_counts_rep"))


fit_dataset <- all_dataset %>%
  data.frame() %>%
  column_to_rownames(var = "GeneID") %>%
  slice(1:5) %>%
  t() %>%
  data.frame %>%
  mutate(group = as.factor(c(rep("0", 20), rep("1", 20))))

### add simulated groups ----------------------------------
probs <- apply(fit_dataset[, -6]*(c(0.4452565, 0.2025354, 0.250560,  0.2531282,  0.2794663 )), 1, sum) - 29.34
probs <- round(inv.logit(mat2),3)
group <- ifelse(probs>0.5, "condition", "control") %>% unlist() %>% data.frame() %>% rownames_to_column(var = "GeneID")
colnames(group) <- c("GeneID", "group")

# save metadata
write_csv(group, file = here("data/simulated_data/metadata_rep"))
