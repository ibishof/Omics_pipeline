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
#install.packages("OOBCurve")
library(OOBCurve)

data <- read_csv("./data/ms_metadata_20190910.csv")
names(data)[names(data) == "Diagnosis_Main"] <- "diagnosis"
names(data)[names(data) == "SMF_ID_or_STS2_Identifier"] <- "sampleid"

data <- data %>% 
  filter(model_cohort=="training") %>% 
  select(sampleid,diagnosis) %>% 
  arrange(sampleid)

# INPUTS POTENTIALLY GOIONG TO COMMAND LINE GO HERE
remove_outliers <- "FALSE"
path1 <- "./scripts/diagnosis/diagnosis_ratios_forest_iter0_"
pull_iter <- "./scripts/diagnosis/iter0"
push_iter <- "./scripts/diagnosis/iter1"
parseCommandArgs()

y <- data$diagnosis
grab_oob <- function(which_mod){
  full_path <- paste(path1,which_mod,path2,sep="")
  mod <- readRDS(full_path)
  #auc_result <- read.delim(full_path, header = FALSE)[1,1]
  return(mod$oobauc)
}

grab_imp <- function(which_mod){
  full_path <- paste(path1,which_mod,path2,sep="")
  mod <- readRDS(full_path)
  return(mod$variable.importance)
}

grab_imp_zero <- function(which_mod){
  full_path <- paste(path1,which_mod,path2,sep="")
  mod <- readRDS(full_path)
  return(mean(mod$variable.importance==0))
}

top_n <- function(x,n=500){
  len<-length(x)
  if(len < n){stop("n must be less than or equal to length of vector.")}
  return(ifelse(x>quantile(x,1-n/len),"yes","no"))
}

pick_top_n <- function(x,n=500){
  return(which(top_n(x,n)=="yes"))
}

path2 <- ".rds"
mod_list <- 1:10
prop_top <- 0.9

split_oob <- sapply(mod_list,grab_oob)
write.table(split_oob,paste(pull_iter,"_oob_results.txt",sep=""),row.names=FALSE)

split_imp <- sapply(mod_list,grab_imp)
avg_imp <- apply(split_imp,1,mean)
num_keep <- ceiling(length(avg_imp)*prop_top)
top_rats <- names(avg_imp)[pick_top_n(avg_imp,num_keep)]
writeLines(top_rats,paste(push_iter,"_ratios.txt",sep=""))
for(i in mod_list){file.remove(paste(path1,i,path2,sep=""))}