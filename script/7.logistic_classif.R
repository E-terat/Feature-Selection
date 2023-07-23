library(stringr)
library(tidyverse)
library(tidymodels)
library(ggplot2)
library(finetune)
library(rsample)
library(parsnip)
library(yardstick)
library(workflows)
library(here)
library(caret)
library(DALEX)

## This .R file is optimized to perform all the logistic models in just one
## 'for' sequence
## It also save all the results in independent folders

analyses_all <- rbind(expand_grid(ds = c("BAS", "ElNet"), d = c("mm"), dea = c("des", "lim"), sub = c(1,0.5,0.2)) %>%
                        data.frame(),
                      expand_grid(ds = c("BAS", "ElNet"), d = c("sim"), dea = c("."),sub = c(1)) %>%
                        data.frame())
                      

for (i in 1:nrow(analyses_all)) {

# FEATURE SELECTION METHOD CAN BE SET TO "BAS" OR "ElNet
feature_selection_method <- analyses_all$ds[i]
# DATASET SELECTION CAN BE SET TO "mm", "hc" or "sim
dataset_selection <- analyses_all$d[i]
# DEA METHOD CAN BE SET TO "des" OR "lim"
dea_method <- ifelse(dataset_selection == "sim", "all", analyses_all$dea[i])
# FEATURE SUBSET (if required to evaluate different subsets -based in feature importance-)
sub <- as.numeric(analyses_all$sub[i])
### metadata file must have samples in rows and info in columns
### norm counts file must have samples in columns and genes in rows


## Metadata loading
metadata_file_d <- ifelse(dataset_selection == "sim",
                          here("data/simulated_data/metadata") , 
                          paste0(here("data/metadata/metadata_"), dataset_selection,"_discovery.csv"))

metadata_file_v <- ifelse(dataset_selection == "sim",
                          here("data/simulated_data/metadata") , 
                          paste0(here("data/metadata/metadata_"), dataset_selection,"_validation.csv"))
metadata_d <- read.csv(metadata_file_d)
metadata_v <- read.csv(metadata_file_v)

colnames(metadata_d) <- tolower(colnames(metadata_d))
colnames(metadata_v) <- tolower(colnames(metadata_v))

## FS results loading
# for simulated data, no dea analysis was performed.
# feature selection methods were directly applyied due to the lower gene number simulated (~ 105)
if (dataset_selection == "sim"){
  fs_file <- paste0(here("results/"),"/", feature_selection_method, "_fit_results/", feature_selection_method, "_" ,dataset_selection, dea_method, ".csv")
  
  fs_results <- ifelse(feature_selection_method == "BAS",
                       read.csv(fs_file) %>% 
                         data.frame() %>%
                         dplyr::slice(-1) %>% 
                         filter(`P.B....0...Y.` > 0.1) %>%
                         dplyr::select(ID),
                       read.csv(fs_file) %>% 
                         data.frame() %>% 
                         filter(Overall > 50) %>%
                         dplyr::select(GeneID)) 
                       
  fs_results <- unlist(fs_results)
  names(fs_results) <- NULL
  fs_results <- fs_results[1:floor(sub*length(fs_results))] # Unmute if a subset of genes is required (for the study 0.5 and 0,2 subsets)
  # real datasets might need a previous dea to pre-select genes in the dataset 
  # so here the data-loading for real data jumps right into the 'else' condition
}else{
  fs_file <- paste0(here("results"),"/", feature_selection_method, "_fit_results/", feature_selection_method, "_", dea_method, "_", dataset_selection, ".csv")
  if (feature_selection_method == "BAS"){                ##gene selection when analyzing bas feature selection method
    fs_results <- read.csv(fs_file) %>% data.frame()   ##just to remove last 5 rows in the file, containing metrics from bas fs
    if ((nrow(fs_results) -6) < nrow(metadata_d)){
      fs_results <- read.csv(dea_file) %>%
        dplyr::select(ID)  %>%
        dplyr::slice(-1) %>% 
        filter(row_number() <= n()-5)%>%
        unlist() 
    }else{
      fs_results <- read.csv(fs_file) %>%
        dplyr::select(ID)  %>%
        dplyr::slice(-1) %>% 
        dplyr::slice(1:(nrow(metadata_d)-2)) %>%
        unlist() 
    }
    names(fs_results) <- NULL
    fs_results <- str_replace(str_replace(str_replace(fs_results,"`", ""),"`", ""), "-", ".")
    
  } else{                                               ##gene selection when analysing ElNet feature selection method
    fs_results <- read.csv(fs_file) %>% data.frame() 
    if ((nrow(fs_results) -5) < nrow(metadata_d)){
      fs_results <- read.csv(fs_file) %>%
        dplyr::select(GeneID)  %>%
        unlist() 
    }else{
      fs_results <- read.csv(fs_file) %>%
        dplyr::select(GeneID)  %>%
        dplyr::slice(1:(nrow(metadata_d)-2)) %>%
        unlist() 
    }
    names(fs_results) <- NULL
    fs_results <- str_replace(str_replace(str_replace(fs_results,"`", ""),"`", ""), "-", ".")
  }
  fs_results <- fs_results[1:floor(sub*length(fs_results))] # Unmute if a subset of genes is required (for the study 0.5 and 0,2 subsets)
}

## (VSN) Normalized counts loading 
norm_counts_file <- paste0(here("data/norm_corrected_data/corr_norm_counts_"), dataset_selection ,"_discovery.csv")
norm_counts <- read_csv(norm_counts_file) %>% 
  column_to_rownames(var = "GeneID")
rownames(norm_counts) <- str_replace(str_replace(str_replace(rownames(norm_counts),"`", ""),"`", ""), "-", ".")
norm_counts_file <- paste0(here("data/norm_corrected_data/corr_norm_counts_"), dataset_selection ,"_validation.csv")
norm_counts_v <- read_csv(norm_counts_file) %>% 
  column_to_rownames(var = "GeneID")
rownames(norm_counts_v) <- str_replace(str_replace(str_replace(rownames(norm_counts_v),"`", ""),"`", ""), "-", ".")

d <- ifelse(dataset_selection == "sim",
            metadata_d %>% filter(sampleid %in% colnames(norm_counts)),
            metadata_d %>% filter(sample %in% colnames(norm_counts)) %>% dplyr::select(group))

d <- d[[1]]

# selection of de genes (dea+fs for real dataset, only fs for simulated dataset)
if (dataset_selection == "sim"){
  discovery <- norm_counts[fs_results,] %>% 
    t() %>% 
    data.frame() %>%
    rownames_to_column(var = "sampleid") %>%
    merge(., metadata_d %>% filter(sampleid %in% colnames(norm_counts)), by = "sampleid") %>%
    mutate(split = "d")
}else{
  discovery <- norm_counts[fs_results,] %>% 
    t() %>% 
    data.frame() %>%
    rownames_to_column(var = "sample") %>%
    merge(., metadata_d %>% dplyr::select(sample, group) %>% filter(sample %in% colnames(norm_counts)), by = "sample") %>%
    mutate(split = "d")
}

v <- ifelse(dataset_selection == "sim",
            metadata_v %>% filter(sampleid %in% colnames(norm_counts_v)) %>% dplyr::select(group) ,
            metadata_v %>% filter(sample %in% colnames(norm_counts_v)) %>% dplyr::select(group))
v <-  v[[1]]
if (dataset_selection == "sim"){
  validation <- norm_counts_v[fs_results,] %>% 
    t() %>% 
    data.frame() %>%
    rownames_to_column(var = "sampleid") %>%
    merge(., metadata_v %>% filter(sampleid %in% colnames(norm_counts_v)), by = "sampleid") %>%
    mutate(split = "v")
  
}else{
  validation <- norm_counts_v[fs_results,] %>% 
    t() %>% 
    data.frame() %>%
    rownames_to_column(var = "sample") %>%
    merge(., metadata_v %>% dplyr::select(sample, group) %>% filter(sample %in% colnames(norm_counts_v)), by = "sample") %>%
    mutate(split = "v")
}


set.seed(123)
model_split <- rsample::group_initial_split(rbind(discovery,validation), split)
disco <- training(model_split) %>% 
  dplyr::select(-split)
valid <- testing(model_split) %>% 
  dplyr::select(-split)



metrics_file_name <- here("results/classifiers/model_metrics_logit.tsv")
metrics_tsv <- file.exists(metrics_file_name)

######### LOGIT n PREDS--------------------------------------------------------------

file_name <- ifelse(dataset_selection == "sim", 
                  paste0("logit_",dataset_selection, dea_method),
                  paste0("logit_",dataset_selection,"_", dea_method,"_", feature_selection_method))


ctrl <- trainControl(method = "LOOCV")  # 10-fold cross-validation
logit_model <- train(group ~ ., data = disco %>% dplyr::select(-1), 
                     method = "glm",
                     family = "binomial", 
                     trControl = ctrl, 
                     metric = "Accuracy")
# library(Boruta)
# boruta_output <- Boruta((as.factor(group)) ~ ., data = data.frame(disco %>% dplyr::select(-1)), doTrace = 0)
# rough_fix_mod <- TentativeRoughFix(boruta_output)
# boruta_signif <- getSelectedAttributes(rough_fix_mod)
# boruta_signif
# importances <- attStats(rough_fix_mod)
# importances <- importances[importances$decision != "Rejected", c("meanImp", "decision")]
# importances[order(-importances$meanImp), ]
# plot(boruta_output, ces.axis = 0.7, las = 2, xlab = "", main = "Feature importance")

sum_model <- summary(logit_model)
predictions <- bind_cols(group = valid$group, 
                         pred.group = predict(logit_model, valid %>% dplyr::select(-1, -group)), 
                         round(predict(logit_model, valid[,-1], type = "prob"),3))

######### METRICS --------------------------------------------------------------

roc_auc <- predictions %>% 
  mutate(group = as.factor(group)) %>%
  yardstick::roc_auc(group, 3)
sens <-  predictions %>% 
  mutate(group = as.factor(group)) %>%
  yardstick::sens(group, pred.group)
spec <- predictions %>% 
  mutate(group = as.factor(group)) %>%
  yardstick::spec(group, pred.group)
ppv <- predictions %>% 
  mutate(group = as.factor(group)) %>%
  yardstick::ppv(group, pred.group)
npv <- predictions %>% 
  mutate(group = as.factor(group)) %>%
  yardstick::sens(group, pred.group)
acc <- predictions %>%
  mutate(group = as.factor(group)) %>%
  accuracy(truth = (group), estimate = pred.group)
roc_curve <- predictions %>%
  mutate(group = as.factor(group)) %>%
  yardstick::roc_curve((group), 3) %>%
  ggplot2::autoplot()

conf_mat <- confusionMatrix(predictions$pred.group, as.factor(valid$group))

aic <- sum_model$aic
dev <- sum_model$deviance
dev.null <- sum_model$null.deviance

#### SAVE METRICS---------------------------------------------------------------
model <- "logit"
test_roc <- roc_auc
test_acc <- acc
test_sens <- sens
test_spec <- spec
test_ppv <- ppv
test_npv <- npv
n_features <- ncol(disco)-2
features <- colnames(disco %>% dplyr::select(-group))[-1]
model_name <- paste(model,
                  feature_selection_method,
                  dataset_selection,
                  dea_method,
                  n_features,
                  sep = "_"
)

metrics_table <- tibble(model_name,
                      model,
                      feature_selection_method,
                      dataset_selection,
                      dea_method,
                      test_roc=test_roc$.estimate,
                      test_acc=test_acc$.estimate,
                      test_spec=test_spec$.estimate,
                      test_sens=test_sens$.estimate,
                      test_ppv=test_ppv$.estimate,
                      test_npv=test_npv$.estimate,
                      aic,
                      dev,
                      dev.null,
                      n_features,
                      features =paste0(features, collapse = "||"))

if (metrics_tsv) {
  metrics_file = read_delim(metrics_file_name,
                          delim = "\t") %>%
    bind_rows(metrics_table) %>%
    write_delim(file = metrics_file_name,
              delim = "\t")
} else {
  if (!dir.exists( here("results/classifiers/logit"))){
    dir.create( here("results/classifiers/logit"))
  }
  write_delim(metrics_table,
            file = metrics_file_name,
            delim = "\t")
}
if (!dir.exists(here(paste0("results/classifiers/logit/", i)))){
  dir.create(here(paste0("results/classifiers/logit/", i)))
}

### PLOTS ---------------------------------------------------------------------
explainer_logit <- DALEX::explain(logit_model,
                                 data = disco %>% dplyr::select(-1) %>% dplyr::select(-group),
                                 y = as.numeric(as.factor(disco$group))-1,
                                 label = "Logistic regression")
diagnostics_plot <- explainer_logit %>% 
  DALEX::model_diagnostics() %>% 
  plot(variable = "y", yvariable = "residuals", smooth = FALSE)

histresiduals_plot <- plot(DALEX::model_performance(explainer_logit), geom="histogram")

featureimp_plot <- explainer_logit %>% model_parts() %>% plot(show_boxplots = FALSE) + 
  ggtitle("Feature Importance ", "") +
  theme(text = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        axis.text = element_text(color = "black")) 

confmat_plot <- cvms::plot_confusion_matrix(as_tibble(conf_mat$table),
                            target_col = "Reference",
                            prediction_col = "Prediction",
                            counts_col = "n")
## SAVE PLOTS------------------------------------------------------------------
ggsave(filename = paste0(model_name, "_residuals.svg"), 
       plot = histresiduals_plot,
       path = here(paste0("results/classifiers/logit/", i)),
       device = "svg", 
       fix_text_size = F)
ggsave(filename = paste0(model_name, "_diagnostics.svg"), 
       plot = diagnostics_plot,
       path = here(paste0("results/classifiers/logit/", i)),
       device = "svg",
       fix_text_size = F)
ggsave(filename = paste0(model_name, "_featureimp.svg"), 
       plot = featureimp_plot,
       path = here(paste0("results/classifiers/logit/", i)),
       device = "svg",
       fix_text_size = F)
ggsave(filename = paste0(model_name, "_confmat.svg"), 
       plot = confmat_plot,
       path = here(paste0("results/classifiers/logit/", i)),
       device = "svg",
       fix_text_size = F)

plot_list <- list()
for (j in 1: nrow(valid)){
  idsample <- (valid %>% dplyr::select(1))[j,]
  nobs <- (valid %>% dplyr::select(-1) %>% dplyr::select(-group))[j,]
  plot_to_save <- list(explainer_logit %>% 
    DALEX::predict_parts(new_observation = nobs) %>% 
    plot() +
    ggtitle(paste0(idsample), "") +
    theme(text = element_text(color = "black"),
          plot.title = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          axis.text = element_text(color = "black")))
  plot_list <- c(plot_list, plot_to_save)
  
  # unmute if individual plots are required
  # ggsave(filename = paste0("predplot_", idsample, ".svg"), 
  #        plot = plot_to_save,
  #        path = here(paste0("results/classifiers/logit/", model_name)),
  #        device = "svg",
  #        fix_text_size = F, 
  #        limitsize = F,
  #        width = 7, height = 5)
  
}
ifelse(nrow(valid) == 18,
       ggplot2::ggsave(filename = "plotpred.svg", 
                plot = ggpubr::ggarrange(plotlist = plot_list, ncol = 3, nrow = 6),
                path = here(paste0("results/classifiers/logit/", i)),
                device = "svg", 
                fix_text_size = F, 
                width = 60, height = 80, units = "cm", 
                bg = "white"),
       ggplot2::ggsave(filename = "plotpred.svg", 
                       plot = ggpubr::ggarrange(plotlist = plot_list, ncol = 2, nrow = 5),
                       path = here(paste0("results/classifiers/logit/", i)),
                       device = "svg", 
                       fix_text_size = F, 
                       width = 40, height = 60, units = "cm", 
                       bg = "white"))

if(dataset_selection == "mm"){ 
  gg <- ggplot(predictions, aes(x=`pred.group`,y=`Multiple Myeloma`, col = group)) +
    geom_jitter(width = 0.2) +  
    theme_minimal() +
    scale_color_manual(labels=c('Mieloma Multiple', 'Control'), values = c("#FFA8F7", "#B0E0F9"))+
    scale_x_discrete(labels=c("Mieloma Multiple", "Control")) +
    xlab("Predicción") +
    ylab("Probabilidad de enfermedad")
  
  ggsave(filename = paste0("pred_all.svg"), 
       plot = gg,
       path = here(paste0("results/classifiers/logit/", i)),
       device = "svg", 
       fix_text_size = F)
}else{
  gg <- ggplot(predictions, aes(x=`pred.group`,y=`condition`, col = group)) +
    geom_jitter(width = 0.2) +  
    theme_minimal() +
    scale_color_manual(labels=c('Enfermedad', 'Control'), values = c("#FFC0A8", "#B0E0E6"))+
    scale_x_discrete(labels=c("Enfermedad", "Control")) +
    xlab("Predicción") +
    ylab("Probabilidad de enfermedad")
  
  ggsave(filename = paste0("pred_all.svg"), 
       plot = gg,
       path = here(paste0("results/classifiers/logit/", i)),
       device = "svg", 
       fix_text_size = F)
}

## SAVE MODEL
saveRDS(final_logit_model,
      file = here(paste0("results/classifiers/logit/", i, "/", model_name, ".rds")))
saveRDS(disco,
      file = here(paste0("results/classifiers/logit/", i, "/data_test.rds")))
saveRDS(valid,
      file = here(paste0("results/classifiers/logit/", i, "/data_train.rds")))

}
