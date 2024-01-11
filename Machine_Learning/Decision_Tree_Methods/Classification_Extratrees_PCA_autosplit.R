# First time running this pipeline you will need to install some R packages
# Extratrees is similar to random froest, it can have better performance on noise weak signal data
# Various R packages that are required
library(readr)
library(dplyr)
library(stringr)
library(ranger)
library(lubridate)
library(pROC)
library(ggfortify)
library(janitor)
library(ggplot2)
library(gridExtra)

# Only the first iteration start here
# Load in data that will be used to train
#Autosplit data
if (isFALSE(exists("random_forest_input"))){
  setwd("C:\\Users\\file\\path")
  data <- read.csv("clean_data.csv", check.names=FALSE)
  data <- data[,-1]
  data_size=nrow(data)
  indexes=sample(1:nrow(data),size = 0.8*data_size)
  training=data[indexes,]
  validation=data[-indexes,]
  
  #This line was added to have the pipeline work with categorical data
  training$Diagnosis = as.factor(training$Diagnosis)
  # Load in validation dataset
  validation$Diagnosis = as.factor(validation$Diagnosis)
}

#####################################################################
##### All iterations after the first start here
#####################################################################

# Model building and feature selection will end when 3 columns aka 2 features are left
while (as.numeric(ncol(training)) >3 ){
  if (exists("random_forest_input")){
    training <- random_forest_input
  }
  # Selecting tuning parameter based on the size of the data
  # if() statement was added to account for iterations with very few ratios
  training <- as.data.frame(training)
  nvars <- floor(3*sqrt(dim(training)[2]))
  if(nvars >= dim(training)[2]-1){nvars <- floor(sqrt(dim(training)[2]))}
  
  
  # Create the model using ranger
  mod <- ranger(data=training,           # The full dataset containing predictors and response variables
                num.trees=min(30000,     # number of trees in the forest
                length(training)),
                splitrule = "extratrees", #If we use splitrule = "extratrees", replace = FALSE, sample.fraction = 1)  we are using Geurts et al. implementation.
                replace = FALSE,
                sample.fraction = 1,
                max.depth = 8,
                min.node.size = 3,
                mtry=nvars,              # number of variables to sample for each split
                importance='impurity',   # Type of variable importance measure to calculate
                seed= 87,                # setting a seed for reproducibility
                dependent.variable.name="Diagnosis",  
                classification = TRUE
                )
  
  # Now that we built our model lets check it's accuracy and save it
  # Check accuracy
  train_predictions <- predict(mod,training)
  auc_train<- roc(as.numeric(training$Diagnosis), as.numeric(as.factor(train_predictions$predictions)))
  
  # Use below for more than two classes
  #auc_train<- roc(as.numeric(training$Diagnosis), as.numeric(as.factor(mod$predictions)))
  # Use if more than one category
  #auc_train<- multiclass.roc(as.numeric(training$Stage), as.numeric(as.factor(mod$predictions)))
  # View results
  auc_train$auc
  mod$confusion.matrix
  
  # Create prediction for validation dataset
  # Calculate accuracy of validation predictions using auc
  predictions <- predict(mod,validation)
  auc_test<- roc(as.numeric(validation$Diagnosis), as.numeric(as.factor(predictions$predictions)))
  auc_test$auc
  
  # For multiple classes
  # test_auc_multi<- multiclass.roc(as.numeric(validation$Diagnosis), as.numeric(as.factor(predictions$predictions)))
  # test_auc_multi$auc
  
  
  # Add everything to list
  if (isFALSE(exists("list_test_auc"))){
    list_test_auc <- rbind(c( ncol(training), auc_test$auc))
  }
  if (exists("list_test_auc")){
    list_test_auc <- rbind(list_test_auc, c(ncol(training), auc_test$auc))
  }
  # Create training auc list
  if (isFALSE(exists("list_train_auc"))){
    list_train_auc <- rbind(c( ncol(training), auc_train$auc))
  }
  if (exists("list_train_auc")){
    list_train_auc <- rbind(list_train_auc, c(ncol(training), auc_train$auc))
  }
  
  # Add predicton error to list
  if (isFALSE(exists("list_prediction_error"))){
    list_prediction_error <- rbind(c( ncol(training), mod$prediction.error))
  }
  if (exists("list_prediction_error")){
    list_prediction_error <- rbind(list_prediction_error, c(ncol(training), mod$prediction.error))
  }
  
  # Save feature importance list
  csv_file_name <- paste0("importance", "_", as.numeric(ncol(training)), ".csv")
  write.csv(mod$variable.importance, csv_file_name)
  
  
  # Make PCA plot
  # Remove zero variance columns
  input<- training
  input <- input %>%
    clean_names()
  input$Diagnosis = as.numeric(as.factor(training$Diagnosis))
  input <- input[ , which(apply(input, 2, var) != 0)]
  # Calculate PCA
  pca_res <- prcomp(input, scale. = TRUE)
  par(mar=c(5,6,4,2))
  title <- paste0("Training ",as.numeric(ncol(training)))
  plot1 <- autoplot(pca_res,
                    data = training,
                    colour = 'Diagnosis',
                    main = title)
  plot1<- plot1 + theme(legend.position = "bottom")
  #print(plot1)

  # Make PCA plot
  # Remove zero variance columns
  input<- validation
  input <- input %>%
    clean_names()
  input$Diagnosis = as.numeric(as.factor(validation$Diagnosis))
  input <- input[ , which(apply(input, 2, var) != 0)]
  # Calculate PCA
  pca_res <- prcomp(input, scale. = TRUE)
  par(mfrow=c(1,1))
  par(mar=c(5,6,4,2))
  title <- paste0("Test ",as.numeric(ncol(training)))
  plot2 <- autoplot(pca_res,
                    data = validation,
                    colour = 'Diagnosis',
                    main = title)
  plot2<- plot2 + theme(legend.position = "bottom")

# Create the PCA plots
  par(mfrow=c(1,2))
  par(mar=c(5,6,4,2))
  grid.arrange(plot1, plot2, nrow = 1)

  
  
  # Recursive feature elimination
  # Take the top 90% of features and create a new list
  var_importance <- as.data.frame(mod$variable.importance)
  ordered <- var_importance[order(-var_importance$`mod$variable.importance`), , drop = FALSE]
  N <- nrow(ordered)
  n <- ceiling(N*.9)
  top90per <- ordered[1:n, ,drop = FALSE]
  names <- unlist(row.names(top90per))
  
  # Take off the one worst feature if less than 10
  if (nrow(ordered) < 10) {
    ordered <- head(ordered, -1)
    names <- unlist(row.names(ordered))
  }
  
  random_forest_input <-training %>% 
    select(Diagnosis,names)
  
  validation <-validation %>% 
    select(Diagnosis,names)
  
  
  # Save model
  rds_file_name <- paste0("Diagnosis", as.numeric(ncol(training)), ".rds")
  saveRDS(mod, file = rds_file_name)
  
  
}

##################################################
##################################################
# End of pipeline everything below is optional
# Write summary of AUC
write.csv(list_train_auc, "list_train_auc.csv") 
write.csv(list_test_auc, "list_test_auc.csv") 

