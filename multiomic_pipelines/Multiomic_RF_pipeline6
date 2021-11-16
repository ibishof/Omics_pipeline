# The pipeline ingest data, builds models, selects top features, and summarizes performance.

## Automate pipeline
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ranger)
library(caret)

################ Functions ################

# This function removes NA and NaN
data.prep <- function(data) {
  # Remove NA or NaN
  data[is.na(data)] = 0
  is.nan.data.frame <- function(x)
    do.call(cbind, lapply(x, is.nan))
  data[is.nan(data)] <- 0
  data <- data
}

# This funcation reduces transcript table down based on top expression features
transcripts.prep <- function(transcripts) {
  # Subset transcripts
  expdata.dir <- paste0(filePath,"expdata")
  setwd(expdata.dir)
  results <- read.csv("results_summarised.csv")
  results <- results[order(-results$mean.test_cor.),]
  # Select best interation
  best_interation <- results[1,]
  best_interation <- best_interation$V1
  best_interation <- paste0(best_interation, "_variable_importance_calc.csv")
  top_features <- read.csv(best_interation)
  # Clean column names
  features_list <- as.vector(top_features$Features)
  feature_list_clean <- gsub("\\..*","", features_list)
  feature_list_clean <- unique(feature_list_clean)
  # Loop though and select columns with shared partial string
  for (i in feature_list_clean){
    if (exists("focus")){
      new <- grep( i, colnames(transcripts))
      focus <- c(focus, new)
    }
    if (isFALSE(exists("focus"))){
      focus <- grep( i, colnames(transcripts))
    }
  }
  # One to get DepMap_ID
  transcripts_focus <- transcripts[,c(1,focus)]
  transcripts_meta <- merge(transcripts_focus, metadata, by = "DepMap_ID")
  transcripts_meta <- select(transcripts_meta, -CCLE_Name)
  data <- transcripts_meta
  
}


################################################################################
# Begining of Data ingest

# Pull in metadata
setwd("/home/ibishof/ccle")
metadata <- read.csv("metadata.csv")

################ Variable ################
Metadata_variables <- c("DepMap_ID", "CCLE_Name", "primary_disease", "AUC_IDB001516")
filePath <- "/home/ibishof/analysis/"
metadata <- select(metadata, Metadata_variables)


# Merge with Metadata
# Pull in data
setwd("/home/ibishof/ccle/ready")
expdata_ready <- readRDS("expdata_ready.rds")
hot_mutations <- readRDS("hot_mutations_ready.rds")
metabolites <- readRDS("metabolites_ready.rds")
proteins <- readRDS("proteins_ready.rds")
transcripts <- readRDS("transcripts_ready.rds")
signatures <- readRDS("signatures_ready.rds")
RRPA <- readRDS("RRPA_ready.rds")
RRBS <- readRDS("RRBS_ready.rds")
CRISPR_Kon <- readRDS("CRISPR_kon_ready.rds")
cn <- readRDS("cn_ready.rds")
gd <- readRDS("gd_ready.rds")
gene_effect <- readRDS("gene_effect_ready.rds")
setwd("/home/ibishof/ccle")
del_mutations <- readRDS("mutations_deleterious.rds")

# Everything
omic_types <- c("expdata", "cn", "gd", "gene_effect", "CRISPR_Kon", "del_mutations", "metabolites", "proteins", "RRBS", "RRPA", "signatures", "transcripts")

# Merge with Metadata
expdata_meta <- merge(expdata_ready, metadata, by = "DepMap_ID")
cn_meta <- merge(cn, metadata, by = "DepMap_ID")
gene_effect_meta <- merge(gene_effect, metadata, by = "DepMap_ID")
del_mutations_meta <- merge(del_mutations, metadata,   by = "DepMap_ID")
CRISPR_Kon_meta <- merge(CRISPR_Kon, metadata,   by = "DepMap_ID")
metabolites_meta <- merge(metabolites, metadata, by = "DepMap_ID")
transcripts_meta <- merge(transcripts, metadata, by = "DepMap_ID")
proteins_meta <- merge(proteins, metadata, by = "CCLE_Name")
signatures_meta <- merge(signatures, metadata, by = "CCLE_Name")
RRPA_meta <- merge(RRPA, metadata, by = "CCLE_Name")
RRBS_meta <- merge(RRBS, metadata, by = "CCLE_Name")
gd_meta <- merge(gd, metadata, by = "DepMap_ID")


# Remove unwanted metadata columns
for (i in omic_types){
  table_name <- paste0(i,"_meta")
  temp_data <- get(table_name)
  ready_data <- select(temp_data, -CCLE_Name)
  title <- paste0(i,"_meta")
  assign(title, ready_data)
}

# Make folder for future results
for (i in omic_types){
  dir.create(paste0(filePath, i))
}




####### RF Pipeline start here #######
for (d in omic_types){
remove(results_summarised)
remove(list_train_cor)
remove(list_test_cor)
remove(df2)
remove(list_test_cor1)
remove(list_train_cor1)
remove(variable_importance)
remove(variable_importance_new_fold)
remove(var_importance)
remove(training)
remove(validation)
remove(v_fold)
remove(skip)

# Only the first iteration start here
# Load in data that will be used to train
if (isFALSE(exists("training"))){
  if (d != "transcripts"){
  table_name <- paste0(d,"_meta")
  data <- get(table_name)
  }

  # Prep transcripts as described in above function "transcripts.prep"
  if (d == "transcripts"){    
    data <- transcripts.prep(transcripts)
  }
  
  # Prep data as described in above function "data.prep"
  data <- data.prep(data)

  # Perform k-fold data subsetting
  k = 10
  data_size=nrow(data)
  data_var <- data
  for (w in 1:k){
    indexes=sample(1:nrow(data_var),size = 1/k*data_size)
    fold=data[indexes,]
    assign(paste0("fold",w), fold)
    data_var=data_var[-indexes,]
  }
  
}


# Select Training and validation set
data.dir <- paste0(filePath,d)
setwd(data.dir)
while (as.numeric(ncol(data)) > 1 ){
  for (i in 1:k){
    v_fold <- paste0("fold",i )
    validation <- get(v_fold)
    
    v_lines <- as.numeric(row.names(validation))
    training <- data[-v_lines,]
    
    
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
                         seed=87,              # setting a seed for reproducibility, DUMZ argument from Rswarm utility
                         #regularization.usedepth = TRUE,
                         #regularization.factor = 5,
                         dependent.variable.name="AUC_IDB001516")
    
    
    
    # Use model to predict and then test to see how accurate those predictions are
    training_predictions <- predict(rating_mod,training)
    training_cor <- cor(training_predictions$predictions,training$AUC_IDB001516, method = "spearman")
    cal_rmse_train = mean((training_predictions$predictions-training$AUC_IDB001516)^2)^0.5
    
    # Create predictions using the validation dataset and calcualte accuracy
    predictions <- predict(rating_mod,validation)
    test_cor <- cor(predictions$predictions,validation$AUC_IDB001516, method = "spearman")
    cal_rmse = mean((predictions$predictions-validation$AUC_IDB001516)^2)^0.5
    
    # Plot predictions vs true vlaues
    par(mfrow=c(1,1))
    Title <- paste("Scatter Plot", ncol(training))
    plot(predictions$predictions, validation$AUC_IDB001516, main= Title,
         xlab="Predictions", ylab="True AUC", cex.lab = 1.5, pch=19)
    
    # Merge all variable importance data across folds
    if (i != 1){
      variable_importance_new_fold <- as.data.frame(rating_mod$variable.importance)
      new_name <- paste0(ncol(training),"_Features_", "Fold",i)
      variable_importance_new_fold <- variable_importance_new_fold %>% 
        rename(
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
    
    # Record perfomance and results    
    # Add everything to list
    if (exists("list_test_cor")){
      list_test_cor1 <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), test_cor, cal_rmse)
      list_test_cor <- rbind(list_test_cor, list_test_cor1)
    }
    if (isFALSE(exists("list_test_cor"))){
      list_test_cor <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), test_cor, cal_rmse)
    }
    # Create training AUC_IDB001516 list
    if (exists("list_train_cor")){
      list_train_cor1 <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), training_cor, cal_rmse_train)
      list_train_cor <- rbind(list_train_cor, list_train_cor1)
    }
    if (isFALSE(exists("list_train_cor"))){
      list_train_cor <- cbind(ncol(training), paste0(ncol(training), "_", "fold", i), training_cor, cal_rmse_train)
    }
    
    
    ################### Interative Feature Selection ######################
    skip <- "no"
    if (i == k){
      # Calculate feature importance across folds mean and SD
      variable_importance$Rmean <- rowMeans(variable_importance[,-1])
      variable_importance$SD <- apply(variable_importance[,-1], 1, sd)
      variable_importance <- variable_importance %>%
        mutate(
          Signal = Rmean/SD
        )      
      write.csv(variable_importance, paste0(ncol(training),"_variable_importance_calc.csv"))
      
      # Feature selection based on colinearity, removed colinear features with low avg importance
      if (ncol(training) < 20000){
        if(isFALSE(exists("df2"))){
          # Get list of feature in descending order
          variable_importance <- as.data.frame(variable_importance)
          variable_importance = variable_importance[order(-variable_importance$`Rmean`), , drop = FALSE]
          features.ordered <- variable_importance$Features
          
          # Reorder column in input table highest importance comes first
          input <- training[features.ordered]
          # Remove zero variance columns
          input <- input[ , which(apply(input, 2, var) != 0)]
          
          # Find correlation and remove co-linear features
          df2 = cor(input, method = "spearman")
          hc = findCorrelation(df2, cutoff=0.80) # put any value as a "cutoff"
          hc = sort(hc)
          reduced_Data = input[,-c(hc)]
          features_clean <- colnames(reduced_Data)
          skip <- "yes"
          
          # Reduce data if primary disease exist
          if (length(features_clean) != 0 & is.null(data$primary_disease) == FALSE) {
            data <- data %>%
              select(AUC_IDB001516, primary_disease, features_clean)
            
            # Write table
            #write.csv(reduced_Data, "reduced_Data.csv")
            for (x in 1:k){
              this_fold <- paste0("fold",x )
              this_fold <- get(this_fold)
              this_fold <- this_fold %>% 
                select(AUC_IDB001516, primary_disease, features_clean)
              assign(paste0("fold",x), this_fold)
            }
          }
          
          # Reduce data if primary disease is NULL
          if (length(features_clean) != 0 & is.null(data$primary_disease) == TRUE) {
            data <- data %>%
              select(AUC_IDB001516, features_clean)
            
            # Write table
            #write.csv(reduced_Data, "reduced_Data.csv")
            for (x in 1:k){
              this_fold <- paste0("fold",x )
              this_fold <- get(this_fold)
              this_fold <- this_fold %>% 
                select(AUC_IDB001516, features_clean)
              assign(paste0("fold",x), this_fold)
            }
          }
          
        }
      }
      
      
      # Remove bottom 10% of features
      if (ncol(data) > 11){
        if (skip == "no"){
        variable_importance <- as.data.frame(variable_importance)
        ordered <- variable_importance[order(-variable_importance$`Rmean`), , drop = FALSE]
        N <- nrow(ordered)
        n <- ceiling(N*.9)
        top90per <- ordered[1:n, ,drop = FALSE]
        names <- unlist(top90per$Features)
        data <- data %>%
          select(AUC_IDB001516,names)
        
        # Create folds with smaller feature list   
        for (y in 1:k){
          this_fold <- paste0("fold",y )
          this_fold <- get(this_fold)
          this_fold <- this_fold %>% 
            select(AUC_IDB001516,names)
          assign(paste0("fold",y), this_fold)
        }
        }
      }
      
      # Take off the one worst feature if less than 10
      if (ncol(data) <= 11) {
        ordered <- head(ordered, -1)
        names <- unlist(ordered$Features)
        data <- data %>%
          select(AUC_IDB001516,names)
        # Create smaller folds   
        for (z in 1:k){
          this_fold <- paste0("fold",z )
          this_fold <- get(this_fold)
          this_fold <- this_fold %>% 
            select(AUC_IDB001516,names)
          assign(paste0("fold",z), this_fold)
        }
        
      }
      # Return to importance based feature reduction
      
    } # end of k == i
    
  } 
  
  
  
  remove(variable_importance)
  
}
# End of RF pipeline


# Save model
# rds_file_name <- paste0("AUC_IDB001516", as.numeric(ncol(training)), ".rds")
# saveRDS(rating_mod, file = rds_file_name)
###################################
###################################

# Summarize model performance and results
if (ncol(data) == 1){
  list_test_cor1 <- as.data.frame(list_test_cor)
  list_test_cor1$test_cor <- as.numeric(as.character(list_test_cor1$V3))
  
  # Group by feature number
  results_summarised <- list_test_cor1 %>%
    group_by(V1) %>%
    summarize(mean(test_cor))
  
  results_summarisedSD <- list_test_cor1 %>%
    group_by(V1) %>%
    summarize(sd(test_cor))
  
  results_summarised <- merge(results_summarised, results_summarisedSD, by = "V1")
  
  results_summarised <- results_summarised %>%
    mutate(
      Signal = `mean(test_cor)`/`sd(test_cor)`
    )
  
  # Save the results
  write.csv(list_train_cor, "list_train_cor.csv")
  write.csv(list_test_cor, "list_test_cor.csv")
  write.csv(results_summarised, "results_summarised.csv")
}
}





