## Kendall rank correlation coefficient 
library(here)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggtext)

list_des <- read.csv("results/DESeq2_results_mm_discovery.csv") %>%
  filter(diffexpressed == "yes") %>%
  pull(gene)

list_lim <- read.csv("results/limma_results_mm_discovery.csv") %>%
  filter(diffexpressed == "yes") %>%
  pull(gene)

listBAS_lim <- read.csv("results/BAS_fit_results/BAS_lim_mm.csv") %>% 
  select(c(ID, `P.B....0...Y.`)) %>% 
  na.exclude() %>%
  filter(!`P.B....0...Y.`==1) %>% 
  pull(ID) 

listBAS_des <- read.csv("results/BAS_fit_results/BAS_des_mm.csv") %>% 
  select(c(ID, `P.B....0...Y.`)) %>% 
  na.exclude() %>%
  filter(!`P.B....0...Y.`==1)  %>% 
  pull(ID)

listElNet_lim <- read.csv("results/ElNet_fit_results/ElNet_lim_mm.csv") %>% 
  select(c(GeneID, Overall)) %>% 
  pull(GeneID) 

listElNet_des <- read.csv("results/ElNet_fit_results/ElNet_des_mm.csv") %>% 
  select(c(GeneID, Overall)) %>%
  pull(GeneID)

list = list(DESeq2 = list_des,
            Limma = list_lim,
            BAS_lim = listBAS_lim,
            BAS_des = listBAS_des,
            ElNet_lim = listElNet_lim, 
            ElNet_des = listElNet_des)
length(list)
plotcorr_kendall <- function(list_of_lists, color = c("slategray")){
  gglist <- list()
  for (i in 1: (length(list_of_lists)-1)){
    for (j in (i+1):length(list_of_lists)){
      listA <- list_of_lists[[i]]
      nameA <- names(list_of_lists[i])
      listB <- list_of_lists[[j]]
      nameB <- names(list_of_lists[j])
      
      if (length(listA)-length(listB) >= 0){
        listA <- listA[1:length(listB)]
      }else{
        listB <- listB[1:length(listA)]
      }
      
      numA <- seq(listA[1:length(listB)]) 
      numB <- match(listA, listB)
      
      difnum <- abs(numA-numB)
      
      k.test <- cor.test(numA,numB, method="kendall")
      k.res <- c(paste0("tau: ",round(k.test$estimate,3)),
                 paste0("p-value: ",round(k.test$p.value,3)))
      data <- data.frame(numA = numA/max(na.omit(c(numA,numB))), 
                         numB = numB/max(na.omit(c(numA,numB))),
                         col = ifelse(difnum<=5, 1,0),
                         label = c(paste(k.res, collapse = "<br>"), 
                                   rep(NA,length(numA)-1)))
      
      gg <- ggplot(data, aes(x=numA, y = numB, col)) + 
        geom_point(alpha = 0.3, 
                   aes(color =  col), 
                   size = 4) + 
        scale_y_continuous(limits = c(0, 1)) +
        scale_x_continuous(limits = c(0, 1)) +
        geom_abline(intercept = 0, slope = 1, color="black", 
                    linetype="dashed", linewidth=1) +
        xlab(nameA) + 
        ylab(nameB) +
        theme_minimal()  +
        theme(legend.position = "none") +
        geom_richtext(
          data = data,
          aes(label = label),
          x = 0.85, 
          y= 0.3,
          size = 5)
      gglist <- c(gglist, list(gg))
    }
  }
  return(gglist)
}

plots <- plotcorr_kendall(list)

ggplot2::ggsave(filename = "plotcorr_kendall.svg", 
       plot = ggpubr::ggarrange(
         plotlist = plots),
       path = here("results/kendall_plots/"),
       device = "svg", 
       fix_text_size = F, 
       width = 40, height = 40, units = "cm")



listElNet_limopt <- read.csv("results/ElNet_fit_results/ElNet_lim_mm.csv") %>% 
  select(c(GeneID, Overall)) %>% 
  pull(GeneID) 

listElNet_desopt <- read.csv("results/ElNet_fit_results/ElNet_des_mm.csv") %>% 
  select(c(GeneID, Overall)) %>%
  pull(GeneID)

listElNet_lim <- read.csv("results/ElNet_fit_results/ElNet_lim_mmround2.csv") %>% 
  select(c(GeneID, Overall)) %>% 
  pull(GeneID) 

listElNet_des <- read.csv("results/ElNet_fit_results/ElNet_des_mmround2.csv") %>% 
  select(c(GeneID, Overall)) %>%
  pull(GeneID)

list = list(limma_ElNet_opt = listElNet_limopt,
            deseq_ElNet_opt = listElNet_desopt,
            limma_ElNet_0.5 = listElNet_lim, 
            des_ElNet_0.5 = listElNet_des)

plots <- plotcorr_kendall(list)

ggplot2::ggsave(filename = "plotcorr_kendall_elnetopt.svg", 
                plot = ggpubr::ggarrange(
                  plotlist = plots),
                path = here("results/kendall_plots/"),
                device = "svg", 
                fix_text_size = F, 
                width = 40, height = 25, units = "cm")


list_sim_des <- read.csv("results/DESeq2_results_sim_discovery.csv") %>% 
  arrange(desc(log2fc)) %>%
  select(c(gene)) %>% 
  pull(gene) 

list_sim_lim <- read.csv("results/limma_results_sim_discovery.csv") %>% 
  arrange(desc(logFC)) %>%
  select(c(gene)) %>%
  pull(gene)

list = list(list_sim_des = list_sim_des[c(1:5)],
            list_sim_lim = list_sim_lim[c(1:5)])

plots <- plotcorr_kendall(list)

ggplot2::ggsave(filename = "plotcorr_kendall_simulationdes.svg", 
                plot = ggpubr::ggarrange(
                  plotlist = plots),
                path = here("results/kendall_plots/"),
                device = "svg", 
                fix_text_size = F, 
                width = 40, height = 25, units = "cm")
