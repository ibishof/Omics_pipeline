# This script uses random forest for intertative feature selection and then feeds those features to Xgboost and Elastic net. Model performance is recorded at each interation and then summarized at the end.
# Various R packages that are required
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ranger)
library(caret)
library(glmnet)
library(xgboost)

#Only first interation starts here
#Bring in training data
#set working directory
#training <- read.csv("train.csv")
#validation <- read.csv("validation.csv")
remove(random_forest_input)
remove(list_train_cor)
remove(list_test_cor)
remove(df2)
remove(list_test_cor1)
remove(list_train_cor1)
remove(variable_importance)
remove(variable_importance_new_fold)
remove(var_importance)
remove(skip)
remove(training)
remove(validation)
remove(v_fold)

# Prevent Masking
select <- dplyr::select
rename <- dplyr::rename
filter <- dplyr::filter


# Prep data
setwd("~/data/brains_proteomics")
data <- readRDS("data_metadata.rds")
data <- filter(data, CDR < 100)

# Create new feature
data <- select(data, -sample, -Group, -Braak, -MMSE)
# Identify the non-numeric columns
non_numeric_columns <- sapply(data, function(x) !is.numeric(x))
# Convert non-numeric columns to numeric
data[, non_numeric_columns] <- lapply(data[, non_numeric_columns], function(x) as.numeric(as.factor(x)))
# data <- select(data, CDR, X7_variable_importance_calc$Features)

# Filter data
# Remove rows with too many NaN
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data[is.nan(data)] <- NA
data[is.na(data)] = 0

# data <- select(data, CDR, X16_variable_importance_calc$Features)
# Now that data is "clean" Let get ready for model building
# This creates "fold" for k-cross fold validation.
# Set K to the number of folds one desires
# Only the first iteration start here
# Load in data that will be used to train
if (isFALSE(exists("random_forest_input"))){
  # Perform k-fold data subsetting
  k = 5
  data_size=nrow(data)
  data_var <- data
  for (w in 1:k){
    indexes=sample(1:nrow(data_var),size = 1/k*data_size)
    fold=data_var[indexes,]
    assign(paste0("fold",w), fold)
    data_var=data_var[-indexes,]
    what_is_i <- w
  }
  
}



#####################################################################

##### All interations after the first start here
##### Use line below after first interation
##### Leave line below commented out for first intration
#####################################################################

#training <- random_forest_input

# Selecting tuning parameter based on the size of the data
# if() statement was added to account for iterations with very few ratios
setwd("~/data/brains_proteomics/models/CDR")
while (as.numeric(ncol(data)) > 1 ){
  for (i in 1:k){
    v_fold <- paste0("fold",i )
    validation <- get(v_fold)
    
    remove_sample <- row.names(validation)
    training <- data[!rownames(data) %in% remove_sample, ] 
    
    
    # Selecting tuning parameter based on the size of the data
    # if() statement was added to account for iterations with very few ratios
    training <- as.data.frame(training)
    nvars <- floor(3*sqrt(dim(training)[2]))
    if(nvars >= dim(training)[2]-1){nvars <- floor(sqrt(dim(training)[2]))}
    
    
    
    # Create the model using ranger
    rating_mod <- ranger(data=training, # The full dataset containing predictors and response variables
                         num.trees=4000, # number of trees in the forest
                         mtry=nvars,      # number of variables to sample for each split
                         importance='impurity',  # Type of variable importance measure to calculate
                         num.threads=16,      # Number of threads for parallelization
                         write.forest = TRUE,
                         seed=i*10,              # setting a seed for reproducibility, DUMZ argument from Rswarm utility
                         #regularization.usedepth = TRUE,
                         #regularization.factor = 5,
                         dependent.variable.name="CDR")
    
    
    
    # Use model to predict and then test to see how accurate those predictions are
    training_predictions <- predict(rating_mod,training)
    training_cor <- cor(training_predictions$predictions,training$CDR, method = "spearman")
    cal_rmse_train = mean((training_predictions$predictions-training$CDR)^2)^0.5
    
    # Create predictions using the validation dataset and calculate accuracy
    predictions <- predict(rating_mod,validation)
    test_cor <- cor(predictions$predictions,validation$CDR, method = "spearman")
    cal_rmse = mean((predictions$predictions-validation$CDR)^2)^0.5
    mae <- mean(abs(predictions$predictions-validation$CDR))
    medae <-median(abs(predictions$predictions-validation$CDR))
    
    # Plot predictions vs true values
    if (i == 1){
      par(mfrow=c(1,5))
    }
    Title <- paste("Scatter Plot", ncol(training))
    plot(predictions$predictions, validation$CDR, main= Title,
         xlab="Predictions", ylab="True CDR", cex.lab = 1.5, pch=19)
    
    
    ############# XgBoost #################################
    # Prepare training and validation data
    train_data <- xgb.DMatrix(data = as.matrix(training[,-which(names(training) == "CDR")]), label = training$CDR)
    val_data <- xgb.DMatrix(data = as.matrix(validation[,-which(names(validation) == "CDR")]), label = validation$CDR)
    
    # Set parameters for the xgboost model
    params <- list(
      objective = "reg:squarederror", # Objective for regression tasks
      eval_metric = "rmse",           # Evaluation metric: root mean squared error
      eta = 0.1,                      # Learning rate
      max_depth = 6,                  # Maximum depth of a tree
      min_child_weight = 1,
      nthread = 2                     # Number of parallel threads used to run xgboost
    )
    
    # Train the xgboost model
    xgb_model <- xgb.train(params, train_data, nrounds = 100, watchlist = list(train = train_data, val = val_data), early_stopping_rounds = 10, verbose = 1)
    
    # Predict CDR for the validation data and training
    train_pred <- predict(xgb_model, train_data)
    val_pred <- predict(xgb_model, val_data)
    
    # Display the predictions
    train_cor_xgb <- cor(train_pred, training$CDR, method = "spearman")
    test_cor_xgb <- cor(val_pred, validation$CDR, method = "spearman")
    cal_rmse_xgb = mean((val_pred-validation$CDR)^2)^0.5
    mae_xgb <- mean(abs(val_pred-validation$CDR))
    medae_xgb <-median(abs(val_pred-validation$CDR))
    
    
    
    ############# Elastic Net #################################
    # data$rf_age <- training_predictions$predictions
    x <- as.matrix(select(training, -CDR))
    
    # Response Variable
    y <- as.factor(training$CDR)
    y2 <- as.factor(validation$CDR)
    
    # Encode the response variable as a binary outcome
    y <- as.numeric(y) - 1
    y2 <- as.numeric(y2) -1
    
    # Make validation matrix
    newx <- as.matrix(select(validation, -CDR))
    
    # Cross-validation
    cvfit <- cv.glmnet(x, y, alpha=.5)
    
    # Make prediction with new model
    predictions <- predict(cvfit, newx = x, s = "lambda.min")
    
    # Plot predictions 
    corr_enet_train <- cor(training$CDR, predictions)
    
    # Plot Validation prediction
    # Make prediction with new model
    predictions <- predict(cvfit, newx = newx, s = "lambda.min")
    
    # Plot predictions 
    corr_enet_test <- cor(validation$CDR, predictions, use = "complete.obs")
    mae_enet <- mean(abs(predictions-validation$CDR))
    medae_enet <-median(abs(predictions-validation$CDR))
   
    
    # Merge all variable importance data
    if (i != 1){
      variable_importance_new_fold <- as.data.frame(rating_mod$variable.importance)
      new_name <- paste0(ncol(training),"_Features_", "Fold",i)
      variable_importance_new_fold <- variable_importance_new_fold %>% 
        dplyr::rename(
          !!new_name := 'rating_mod$variable.importance')
      
      variable_importance_new_fold$Features <- row.names(variable_importance_new_fold)
      variable_importance <- merge(variable_importance, variable_importance_new_fold, by="Features", all=TRUE)
    }
    
    # For first fold
    if (i == 1){
      variable_importance <- as.data.frame(rating_mod$variable.importance)
      variable_importance$Features <- row.names(variable_importance)
      variable_importance <- variable_importance %>% 
        rename(
          Fold1 = 'rating_mod$variable.importance')
    }
    
    # Control correlation importance based feature removal
    skip <- "no"
    
    # Add everything to list
    if (exists("list_test_cor")){
      list_test_cor1 <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), test_cor, cal_rmse)
      list_test_cor <- rbind(list_test_cor, list_test_cor1)
    }
    if (isFALSE(exists("list_test_cor"))){
      list_test_cor <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), test_cor, cal_rmse)
    }
    # Create training CDR list
    if (exists("list_train_cor")){
      list_train_cor1 <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), training_cor, cal_rmse_train)
      list_train_cor <- rbind(list_train_cor, list_train_cor1)
    }
    if (isFALSE(exists("list_train_cor"))){
      list_train_cor <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), training_cor, cal_rmse_train)
    }
    if (exists("list_test_mae")){
      list_test_mae1 <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), mae)
      list_test_mae <- rbind(list_test_mae, list_test_mae1)
    }
    if (isFALSE(exists("list_test_mae"))){
      list_test_mae <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), mae)
    }
    if (exists("list_test_medae")){
      list_test_medae1 <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), medae)
      list_test_medae <- rbind(list_test_medae, list_test_medae1)
    }
    if (isFALSE(exists("list_test_medae"))){
      list_test_medae <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), medae)
    }
    
    
    # Xgboost put everything in list
    if (exists("list_test_cor_xgb")){
      list_test_cor1_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), test_cor_xgb)
      list_test_cor_xgb <- rbind(list_test_cor_xgb, list_test_cor1_xgb)
    }
    if (isFALSE(exists("list_test_cor_xgb"))){
      list_test_cor_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), test_cor_xgb)
    }
    # Create training CDR list
    if (exists("list_train_cor_xgb")){
      list_train_cor1_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), train_cor_xgb)
      list_train_cor_xgb <- rbind(list_train_cor_xgb, list_train_cor1_xgb)
    }
    if (isFALSE(exists("list_train_cor_xgb"))){
      list_train_cor_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), train_cor_xgb)
    }
    if (exists("list_test_mae_xgb")){
      list_test_mae1_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), mae_xgb)
      list_test_mae_xgb <- rbind(list_test_mae_xgb, list_test_mae1_xgb)
    }
    if (isFALSE(exists("list_test_mae_xgb"))){
      list_test_mae_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), mae_xgb)
    }
    if (exists("list_test_medae_xgb")){
      list_test_medae1_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), medae_xgb)
      list_test_medae_xgb <- rbind(list_test_medae_xgb, list_test_medae1_xgb)
    }
    if (isFALSE(exists("list_test_medae_xgb"))){
      list_test_medae_xgb <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), medae_xgb)
    }
    
    
    # E-net version of everything
    # Add everything to list
    if (exists("list_test_cor_enet")){
      list_test_cor1_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), corr_enet_test)
      list_test_cor_enet <- rbind(list_test_cor_enet, list_test_cor1_enet)
    }
    if (isFALSE(exists("list_test_cor_enet"))){
      list_test_cor_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), corr_enet_test)
    }
    # Create training CDR list
    if (exists("list_train_cor_enet")){
      list_train_cor1_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), corr_enet_train)
      list_train_cor_enet <- rbind(list_train_cor_enet, list_train_cor1_enet)
    }
    if (isFALSE(exists("list_train_cor_enet"))){
      list_train_cor_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), corr_enet_train)
    }
    if (exists("list_test_mae_enet")){
      list_test_mae1_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), mae_enet)
      list_test_mae_enet <- rbind(list_test_mae_enet, list_test_mae1_enet)
    }
    if (isFALSE(exists("list_test_mae_enet"))){
      list_test_mae_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), mae_enet)
    }
    if (exists("list_test_medae_enet")){
      list_test_medae1_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), medae_enet)
      list_test_medae_enet <- rbind(list_test_medae_enet, list_test_medae1_enet)
    }
    if (isFALSE(exists("list_test_medae_enet"))){
      list_test_medae_enet <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), medae_enet)
    }
    
    
    ################### Feature Selection ######################
    if (i == k){
      # Calculate mean and SD
      variable_importance$Rmean <- rowMeans(variable_importance[,-1])
      variable_importance$SD <- apply(variable_importance[,-1], 1, sd)
      variable_importance <- variable_importance %>%
        mutate(
          Signal = Rmean/SD
        )
      
      write.csv(variable_importance, paste0(ncol(training),"_variable_importance_calc.csv"))
      
      # Colinearity filter
      # if (ncol(training) < 20000){
      #   if(isFALSE(exists("df2"))){
      #     # Get list of feature in descending order
      #     variable_importance <- as.data.frame(variable_importance)
      #     variable_importance = variable_importance[order(-variable_importance$`Rmean`), , drop = FALSE]
      #     features.ordered <- variable_importance$Features
      # 
      #     # Reorder column in input table highest importance comes first
      #     input <- training[features.ordered]
      #     # Remove zero variance columns
      #     input <- input[ , which(apply(input, 2, var) != 0)]
      # 
      #     # Find correlation and remove co-linear features
      #     df2 = cor(input, method = "spearman")
      #     hc = findCorrelation(df2, cutoff=0.89) # put any value as a "cutoff"
      #     hc = sort(hc)
      #     reduced_Data = input[,-c(hc)]
      #     features_clean <- colnames(reduced_Data)
      #     skip <- "yes"
      # 
      #     # Reduce data
      #     data <- data %>%
      #       select(CDR,features_clean)
      # 
      #     # Write table
      #     #write.csv(reduced_Data, "reduced_Data.csv")
      #     for (x in 1:k){
      #       this_fold <- paste0("fold",x )
      #       this_fold <- get(this_fold)
      #       this_fold <- this_fold %>%
      #         select(CDR,features_clean)
      #       assign(paste0("fold",x), this_fold)
      #     }
      # 
      #   }
      # }
      
      
      # Remove bottom 10% of features
      if (skip == "no"){
        variable_importance <- as.data.frame(variable_importance)
        ordered <- variable_importance[order(-variable_importance$`Rmean`), , drop = FALSE]
        N <- nrow(ordered)
        n <- ceiling(N*.9)
        top90per <- ordered[1:n, ,drop = FALSE]
        names <- unlist(top90per$Features)
        data <- data %>%
          select(CDR,names)
        
        # Create folds with smaller feature list   
        for (y in 1:k){
          this_fold <- paste0("fold",y )
          this_fold <- get(this_fold)
          this_fold <- this_fold %>% 
            select(CDR,names)
          assign(paste0("fold",y), this_fold)
        }
      }
      # Take off the one worst feature if less than 10
      if (ncol(data) < 11) {
        ordered <- head(ordered, -1)
        names <- unlist(ordered$Features)
        data <- data %>%
          select(CDR,names)
        # Create smaller folds   
        for (z in 1:k){
          this_fold <- paste0("fold",z )
          this_fold <- get(this_fold)
          this_fold <- this_fold %>% 
            select(CDR,names)
          assign(paste0("fold",z), this_fold)
        }
        
      }
      # Return to importance based feature reduction
      skip <- "no"
    } # end of k == i
    
  } 
  
  
  
  remove(variable_importance)
  
}




# Save model
# rds_file_name <- paste0("CDR", as.numeric(ncol(training)), ".rds")
# saveRDS(rating_mod, file = rds_file_name)



###################################
###################################
# End of pipeline everything below is optional
# Optional View variable list
# Save predictions and feature importance

list_test_cor1 <- as.data.frame(list_test_cor)
list_test_cor1$test_cor <- as.numeric(as.character(list_test_cor1$test_cor))

list_test_mae1 <- as.data.frame(list_test_mae)
list_test_mae1$mae <- as.numeric(as.character(list_test_mae1$mae))

list_test_medae1 <- as.data.frame(list_test_medae)
list_test_medae1$medae <- as.numeric(as.character(list_test_medae1$medae))


# XgBoost
list_test_cor1_xgb <- as.data.frame(list_test_cor_xgb)
list_test_cor1_xgb$test_cor_xgb <- as.numeric(as.character(list_test_cor1_xgb$test_cor_xgb))

list_test_mae1_xgb <- as.data.frame(list_test_mae_xgb)
list_test_mae1_xgb$mae_xgb <- as.numeric(as.character(list_test_mae1_xgb$mae_xgb))

list_test_medae1_xgb <- as.data.frame(list_test_medae_xgb)
list_test_medae1_xgb$medae_xgb <- as.numeric(as.character(list_test_medae1_xgb$medae_xgb))

# E-net
list_test_cor1_enet <- as.data.frame(list_test_cor_enet)
list_test_cor1_enet$test_cor_enet <- as.numeric(as.character(list_test_cor1_enet$lambda.min))

list_test_mae1_enet <- as.data.frame(list_test_mae_enet)
list_test_mae1_enet$mae_enet <- as.numeric(as.character(list_test_mae1_enet$mae_enet))

list_test_medae1_enet <- as.data.frame(list_test_medae_enet)
list_test_medae1_enet$medae_enet <- as.numeric(as.character(list_test_medae1_enet$medae_enet))



# Group by feature number
results_summarised <- list_test_cor1 %>%
  group_by(V1) %>%
  summarize(mean(test_cor))

results_summarisedSD <- list_test_cor1 %>%
  group_by(V1) %>%
  summarize(sd(test_cor))

results_summarised_mae <- list_test_mae1 %>%
  group_by(V1) %>%
  summarize(mean(mae))

results_summarised_medae <- list_test_medae1 %>%
  group_by(V1) %>%
  summarize(mean(medae))


results_summarised_xgb <- list_test_cor1_xgb %>%
  group_by(V1) %>%
  summarize(mean(test_cor_xgb))

results_summarised_xgbSD <- list_test_cor1_xgb %>%
  group_by(V1) %>%
  summarize(sd(test_cor_xgb))

results_summarised_mae_xgb <- list_test_mae1_xgb %>%
  group_by(V1) %>%
  summarize(mean(mae_xgb))

results_summarised_medae_xgb <- list_test_medae1_xgb %>%
  group_by(V1) %>%
  summarize(mean(medae_xgb))


results_summarised_enet <- list_test_cor1_enet %>%
  group_by(V1) %>%
  summarize(mean(test_cor_enet))

results_summarised_enetSD <- list_test_cor1_enet %>%
  group_by(V1) %>%
  summarize(sd(test_cor_enet))

results_summarised_mae_enet <- list_test_mae1_enet %>%
  group_by(V1) %>%
  summarize(mean(mae_enet))

results_summarised_medae_enet <- list_test_medae1_enet %>%
  group_by(V1) %>%
  summarize(mean(medae_enet))


results_summarised <- merge(results_summarised, results_summarisedSD, by = "V1")

results_summarised <- merge(results_summarised, results_summarised_mae, by = "V1")

results_summarised <- merge(results_summarised, results_summarised_medae, by = "V1")

results_summarised <- results_summarised %>%
  mutate(
    Signal = `mean(test_cor)`/`sd(test_cor)`
  )



results_summarised <- merge(results_summarised, results_summarised_xgb, by = "V1")

results_summarised <- merge(results_summarised, results_summarised_mae_xgb, by = "V1")

results_summarised <- merge(results_summarised, results_summarised_medae_xgb, by = "V1")

results_summarised$Signal_xgb <- results_summarised$`mean(test_cor_xgb)`/results_summarised_xgbSD$`sd(test_cor_xgb)`



results_summarised <- merge(results_summarised, results_summarised_enet, by = "V1")

results_summarised <- merge(results_summarised, results_summarised_mae_enet, by = "V1")

results_summarised <- merge(results_summarised, results_summarised_medae_enet, by = "V1")

results_summarised$Signal_enet <- results_summarised$`mean(test_cor_enet)`/results_summarised_enetSD$`sd(test_cor_enet)`

# Save the results
write.csv(list_train_cor, "list_train_cor.csv")
write.csv(list_test_cor, "list_test_cor.csv")
write.csv(results_summarised, "results_summarised.csv")


# setwd("~/Rmarkdown/test")
# performance_table <- read.csv("results_summarised.csv")
par(mfrow=c(1,1))
performance_table <- results_summarised
performance_table$V1 <- as.numeric(performance_table$V1)
performance_table$`mean(test_cor)` <- as.numeric(performance_table$`mean(test_cor)`)
performance_table$`sd(test_cor)` <- as.numeric(performance_table$`sd(test_cor)`)
performance_table$Signal <- as.numeric(performance_table$Signal)
avg <- performance_table$`mean(test_cor)`
sdev <- performance_table$`sd(test_cor)`
x <- performance_table$V1

plot(x, avg,
     ylim=range(c(avg-sdev, avg+sdev)),
     xlim=rev(range(c(x))),
     pch=19, xlab="# of Features", ylab="Avg Corr +/- SD",
     main="Performance and Featue Elimination"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)
first_interaction <- filter(performance_table, V1 == max(performance_table$V1))
height <- first_interaction$mean.test_cor.
abline(h=first_interaction$`mean(test_cor)`,col="red", lwd=3, lty=2)

# Make plot for model signal
signal <- performance_table$Signal
plot(x, signal, pch =19, xlab = "Features", ylab = "Avg Cor/SDev", xlim=rev(range(c(x))), main = "Signal and Feature Elmination")
height_signal <- first_interaction$Signal
abline(h=height_signal,col="red", lwd=3, lty=2)



