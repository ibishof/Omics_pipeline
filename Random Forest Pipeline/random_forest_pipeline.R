#Load in libraries
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
source("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\Bibi\\pipeline_for_auroc\\pipeline_for_meeting\\scripts\\pipeline_results_function.R")
setwd('W://')
# Load in data set
stringsAsFactors = FALSE
soma_untreated <- read.csv(".\\data\\random_forest_input.csv",header = TRUE, stringsAsFactors = FALSE)
# Changes column names to better suit pipeline
names(soma_untreated)[names(soma_untreated) == "Diagnosis_Main"] <- "diagnosis"
names(soma_untreated)[names(soma_untreated) == "SMF_ID_or_STS2_Identifier"] <- "sampleid"
soma_untreated$diagnosis[soma_untreated$diagnosis != "MS"]<- "other"

# Ratio fucntion used to select and create ratios for final model
ratio_append <- function(data,ratio_list){
  splitlist <- strsplit(ratio_list,"[_]")
  v1 <- sapply(splitlist,first)
  v2 <- sapply(splitlist,last)
  p <- length(ratio_list)
  for(i in 1:p){
    if(length(splitlist[[i]]) == 2){
      data[,ratio_list[i]] <- data[[v1[i]]]/data[[v2[i]]]
    }else{
      data[,ratio_list[i]] <- data[[v1[i]]]
    }
  }
  return(data)
}

#Selects iteration with the best oob value
diagnosis_plot <- pipeline_results("diagnosis","Diagnosis Baseline Visit",iters=0:43)
diagnosis_plot <- pipeline_results("diagnosis","Diagnosis Baseline Visit",iters=0:100,
                                   best_iter=75)
diagnosis_plot

# Using eye ball
diagnosis_iter <- 40

# Selects features from specified iteration
diagnosis_var <- readLines(paste('./scripts/diagnosis/iter',diagnosis_iter,"_ratios.txt",sep=""))
soma_untreated <- soma_untreated %>% ratio_append(diagnosis_var)

# Creates training and test data set
soma_train <- soma_untreated %>% filter(model_cohort=="training")
soma_test <- soma_untreated %>% filter(model_cohort=="validation")

# Creates nvars and data to be feed into 
nvars <- floor(3*sqrt(length(temp_data)))
temp_data <- as.data.frame(soma_train[,c("diagnosis",diagnosis_var)])
temp_data$diagnosis <- as.factor(temp_data$diagnosis)
diagnosis_mod <- ranger(data=temp_data, num.trees=min(40000,length(temp_data)),mtry=nvars,importance='impurity',
                        seed=225,dependent.variable.name="diagnosis",classification = TRUE)
diagnosis_mod$confusion.matrix

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
msdss_imp_dat <- tibble(ratio = names(diagnosis_mod$variable.importance),
                        importance = as.numeric(diagnosis_mod$variable.importance))
write.csv(msdss_imp_dat, "feature_importance.csv")
