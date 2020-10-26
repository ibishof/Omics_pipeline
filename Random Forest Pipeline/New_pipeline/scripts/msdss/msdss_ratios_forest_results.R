# Various R packages that are requried
library(readr)
library(dplyr)
library(tidyr)
library(bindrcpp)
library(data.table)
library(stringr)
library(ranger)
library(lubridate)
library(parallel)
library(doParallel)
library(batch)

# Load the metadata into an R object that is called "data"
data <- read_csv("./data/ms_metadata_20190910.csv")

# In the following code, we take the metadata, filter down to the training cohort, select only the
# Identifier and the trait of interest (msdss) and convert to a data.table to merge with the somamer ratios. 
# If we wanted to build a model predicting a different outcome with a different subset of data, those steps can
# be included here. GOOD TO MATCH WHAT IS ALREADY SPECIFIED IN MSDSS_RATIOS_FOREST.r SCRIPT
data <- data %>% 
  filter(model_cohort=="training") %>% 
  select(sampleid,msdss) %>% 
  arrange(sampleid)

# Command line input using the batch R package. Tells R where to look for the forest we are analyzing, the 
# ratios that were used in the current iteration, and what to save the next set of ratios as.
# 
path1 <- "./scripts/msdss/msdss_ratios_forest_iter0_"
pull_iter <- "./scripts/msdss/iter0"
push_iter <- "./scripts/msdss/iter1"
parseCommandArgs()

# Create a vector of the outcome/response of interest. Can adapt for different problems.
y <- data$msdss

# Function to grab OOB error from a specific model
grab_oob <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(mod$prediction.error)
}

# Function to grab OOB correlation from a specified model
grab_cor <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(cor(y,mod$predictions))
}

# Function I wrote to calculate the Concordance Correlation Coefficient (CCC) between 
# two variables
ccc_function <- function(x1,x2){
	n <- length(x1)
	xb <- mean(x1)
	yb <- mean(x2)
	sx <- mean((x1-xb)^2)
	sy <- mean((x2-yb)^2)
	sxy <- mean((x1-xb)*(x2-yb))
	out <- (2*sxy)/(sx + sy + (xb - yb)^2)
	return(out)
}

# Function to grab OOB CCC from a specified model
grab_ccc <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(ccc_function(y,mod$predictions))
}

# Function to grab the variable importance measures from a specified model
grab_imp <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(mod$variable.importance)
}

# Function to grab the percentage of variable importance measures that 
# are equal to 0 from a specified model
grab_imp_zero <- function(which_mod){
	full_path <- paste(path1,which_mod,path2,sep="")
	mod <- readRDS(full_path)
	return(mean(mod$variable.importance==0))
}

# Function that I wrote to determine if an element is in the "n" largest elements of the vector.
# Used to pick the top percentage 
top_n <- function(x,n=500){
	 len<-length(x)
	 if(len < n){stop("n must be less than or equal to length of vector.")}
	 return(ifelse(x>quantile(x,1-n/len),"yes","no"))
}

# Using the function above, returns which element of a vector are in the top "n" largest elements
pick_top_n <- function(x,n=500){
	 return(which(top_n(x,n)=="yes"))
}

# Used to grab all the interesting peices in the functions that were written above.
path2 <- ".rds"
mod_list <- 1:10
prop_top <- 0.9

# Grab and save the OOB error from all of the models
split_oob <- sapply(mod_list,grab_oob)
write.table(split_oob,paste(pull_iter,"_oob_results.txt",sep=""),row.names=FALSE)

# Grab and save the OOB correlation from all of the models
split_cor <- sapply(mod_list,grab_cor)
write.table(split_cor,paste(pull_iter,"_oob_cor_results.txt",sep=""),row.names=FALSE)

# Grab and save the OOB concordance from all of the models
split_ccc <- sapply(mod_list,grab_ccc)
write.table(split_ccc,paste(pull_iter,"_oob_ccc_results.txt",sep=""),row.names=FALSE)

# Average together variable importance measures from all models, and return/save those that are in the top 90%
split_imp <- sapply(mod_list,grab_imp)
avg_imp <- apply(split_imp,1,mean)
num_keep <- ceiling(length(avg_imp)*prop_top)
top_rats <- names(avg_imp)[pick_top_n(avg_imp,num_keep)]
writeLines(top_rats,paste(push_iter,"_ratios.txt",sep=""))

# Delete the models so they are not taking up space
for(i in mod_list){file.remove(paste(path1,i,path2,sep=""))}
