## knee plot
library(readr)
library(tidyverse)
library(changepoint)


BAS_simall_final <- read_csv("results/BAS_fit_results/BAS_simall_final.csv") %>%
  arrange(desc(mean)) %>%
  rownames_to_column(var="num") %>%
  mutate(mean = as.numeric(mean), 
         num = as.numeric(num))

BAS_lim_mm_final <- read_csv("results/BAS_fit_results/BAS_lim_mm_final.csv") %>%
  arrange(desc(mean)) %>%
  rownames_to_column(var="num") %>%
  mutate(mean = as.numeric(mean), 
         num = as.numeric(num))

BAS_mm_des_final <- read_csv("results/BAS_fit_results/BAS_mm_des_final.csv") %>%
  arrange(desc(mean)) %>%
  rownames_to_column(var="num") %>%
  mutate(mean = as.numeric(mean), 
         num = as.numeric(num))

cpt <- cpt.var(BAS_simall_final$mean, method = "BinSeg")
cpt@cpts[1]
ggplot(BAS_simall_final, aes(num, mean)) +
  geom_line(color = "gray", size = 0.5) +
  geom_point(color = "steelblue", size = 2, shape = 21, fill = "steelblue") +
  theme_minimal() +
  labs(x = "Variables", y = "Probability", title = "Elbow plot: BAS probs in Simulated data") +
  geom_vline(xintercept = cpt@cpts[1], color = "black", linetype = 2)
ggsave(filename = "BAS_simall_final_elbow.png", device = "png", bg = "white")

cpt <- cpt.var(BAS_lim_mm_final$mean, method = "BinSeg")
cpt@cpts[1]
ggplot(BAS_lim_mm_final, aes(num, mean)) +
  geom_line(color = "gray", size = 0.5) +
  geom_point(color = "steelblue", size = 2, shape = 21, fill = "steelblue") +
  theme_minimal() +
  labs(x = "Variables", y = "Probability", title = "Elbow plot: Limma + BAS probs in Multiple Myeloma data") +
  geom_vline(xintercept = cpt@cpts[1], color = "black", linetype = 2)
ggsave(filename = "BAS_lim_mm_final_elbow.png", device = "png", bg = "white")


cpt <- cpt.var(BAS_mm_des_final$mean, method = "BinSeg")
cpt@cpts[1]
ggplot(BAS_mm_des_final, aes(num, mean)) +
  geom_line(color = "gray", size = 0.5) +
  geom_point(color = "steelblue", size = 2, shape = 21, fill = "steelblue") +
  theme_minimal() +
  labs(x = "Variables", y = "Probability", title = "Elbow plot: DESeq2+BAS probs in Multiple Myeloma data") +
  geom_vline(xintercept = cpt@cpts[1], color = "black", linetype = 2)
ggsave(filename = "BAS_mm_des_final_elbow.png", device = "png", bg = "white")

