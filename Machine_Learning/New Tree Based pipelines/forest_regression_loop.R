# Various R packages that are requried
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ranger)


#Only first interation starts here
#Bring in training data
#set working directory
training <- read.csv("train.csv")
validation <- read.csv("validation.csv")
training <- training_cons
validation <- validation_cons
remove(random_forest_input)
# Only the first iteration start here
# Load in data that will be used to train
if (isFALSE(exists("random_forest_input"))){
  setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\metabolite_regression\\ratios")
  data <- read.csv("inputz.csv", check.names=FALSE)
  data <- data[,-1]
  data_size=nrow(data)
  indexes=sample(1:nrow(data),size = 0.75*data_size)
  training=data[indexes,]
  validation=data[-indexes,]
  # For later use
  training_cons <- training
  validation_cons <- validation

}

#####################################################################
##### All iterations after the first start here
#####################################################################

# Selecting tuning parameter based on the size of the data
# if() statement was added to account for iterations with very few ratios
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
rating_mod <- ranger(data=training, # The full dataset containing predictors and response variables
              num.trees=40000, # number of trees in the forest
              mtry=nvars,      # number of variables to sample for each split
              importance='impurity',  # Type of variable importance measure to calculate
              num.threads=2,      # Number of threads for parallelization
              #write.forest = TRUE,
	            seed=10943905,              # setting a seed for reproducibility, DUMZ argument from Rswarm utility
              splitrule = "extratrees",
              max.depth = 8, # Maximal tree depth. A value of NULL or 0 (the default) corresponds to unlimited depth, 1 to tree stumps (1 split per tree).
              min.node.size = 3,
              sample.fraction = 1,
              replace = FALSE,
              #regularization.usedepth = TRUE,
              #regularization.factor = 5,
              keep.inbag = TRUE,
	            dependent.variable.name="Diagnosis")  # Telling ranger what the response variable is



# Use model to predict and then test to see how accurate those predictions are
training_predictions <- predict(rating_mod,training)
training_cor <- cor(training_predictions$predictions,training$Diagnosis, method = "pearson")


# Create predictions using the validation dataset and calcualte accuracy
predictions <- predict(rating_mod,validation)
test_cor <- cor(predictions$predictions,validation$Diagnosis, method = "pearson")

# Plot predictions vs true vlaues
plot(predictions$predictions, validation$Diagnosis, main="Scatter Plot",
     xlab="Predictions", ylab="True Stage", cex.lab = 1.5, pch=19)


# Add everything to list
if (isFALSE(exists("list_test_cor"))){
  list_test_cor <- rbind(c( ncol(training), test_cor))
}
if (exists("list_test_cor")){
  list_test_cor <- rbind(list_test_cor, c(ncol(training), test_cor))
}
# Create training auc list
if (isFALSE(exists("list_train_cor"))){
  list_train_cor <- rbind(c( ncol(training), training_cor))
}
if (exists("list_train_cor")){
  list_train_cor <- rbind(list_train_cor, c(ncol(training), training_cor))
}


# Save feature importance list
csv_file_name <- paste0("importance", "_", as.numeric(ncol(training)), ".csv")
write.csv(rating_mod$variable.importance, csv_file_name)



# Take the top 90% of features and create a new list
var_importance <- as.data.frame(rating_mod$variable.importance)
ordered <- var_importance[order(-var_importance$`rating_mod$variable.importance`), , drop = FALSE]
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
rds_file_name <- paste0("Regression", as.numeric(ncol(training)), ".rds")
saveRDS(rating_mod, file = rds_file_name)


}

###################################
###################################
# End of pipeline everything below is optional
# Optional View variable list
# Save predictions and feature importance
View(rating_mod$variable.importance)
write.csv(list_train_cor, "list_train_cor.csv")
write.csv(list_train_cor, "list_train_cor.csv")



