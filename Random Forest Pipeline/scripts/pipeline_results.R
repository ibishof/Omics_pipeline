library(readr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(bindrcpp, quietly = TRUE)
library(stringr, quietly = TRUE)
library(ranger)
library(ggplot2)

# Load function that does the heavy lifting
source("./scripts/pipeline_results_function.R")

soma_untreated <- read_csv("./data/ms_metadata_20190910.csv")

names(soma_untreated)[names(soma_untreated) == "Diagnosis_Main"] <- "diagnosis"
names(soma_untreated)[names(soma_untreated) == "SMF_ID_or_STS2_Identifier"] <- "sampleid"
soma_untreated$diagnosis[soma_untreated$diagnosis != "MS"] <- "Non-MS"

ratio_append <- function(data,ratio_list){
  splitlist <- strsplit(ratio_list,"[_]")
  v1 <- sapply(splitlist,first)
  v2 <- sapply(splitlist,last)
  p <- length(ratio_list)
  for(i in 1:p){
    data[,ratio_list[i]] <- data[[v1[i]]]/data[[v2[i]]]
  }
  return(data)
}

diagnosis_plot <- pipeline_results("diagnosis","Diagnosis Baseline Visit",iters=0:104)
diagnosis_plot <- pipeline_results("diagnosis","Diagnosis Baseline Visit",iters=0:90,
                                   best_iter=71)
diagnosis_plot

# Using eye ball
diagnosis_iter <- 71


diagnosis_var <- readLines(paste('./scripts/diagnosis/iter',diagnosis_iter,"_ratios.txt",sep=""))
soma_untreated <- soma_untreated %>% ratio_append(diagnosis_var)

soma_train <- soma_untreated %>% filter(model_cohort=="training")
soma_test <- soma_untreated %>% filter(model_cohort=="validation")

nvars <- floor(3*sqrt(length(diagnosis_var)))
temp_data <- as.data.frame(soma_train[,c("diagnosis",diagnosis_var)])
temp_data$diagnosis <- as.factor(temp_data$diagnosis)
diagnosis_mod <- ranger(data=temp_data, num.trees=min(40000,length(diagnosis_var)),mtry=nvars,importance='impurity',
                        seed=223,dependent.variable.name="diagnosis",classification = TRUE)

######################################################################
##################   Add predictions to Dataset ######################
######################################################################
soma_train$soma_diagnosis <- diagnosis_mod$predictions
soma_test$soma_diagnosis <- predict(diagnosis_mod,soma_test)$predictions

# AUROC
library(pROC)
roc(as.numeric(soma_test$soma_diagnosis), as.numeric(as.factor(soma_test$diagnosis)))

# Confusion Matrix
soma_untreated <- bind_rows(soma_train,soma_test)
table(soma_untreated$diagnosis, soma_untreated$soma_diagnosis)

############################################################################
############### Ratios and Variable Importance in Models ###################
############################################################################

# Create a dataset that contains the ratios that were used in the final model
# along with the variable importance measures for said ratios
msdss_imp_dat <- tibble(ratio = names(msdss_mod$variable.importance),
                        importance = as.numeric(diagnosis_mod$variable.importance))
write.csv(msdss_imp_dat, "msdss_imp_dat.csv")
