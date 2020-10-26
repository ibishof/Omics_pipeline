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

# Biowulf function to determine the correct number of CPUS to use
detectBatchCPUs <- function() { 
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
  if (is.na(ncores)) { 
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
  } 
  if (is.na(ncores)) { 
    return(4) # for helix
  } 
  return(ncores) 
}

# Load the metadata into an R object that is called "data"
data <- read_csv("./data/ms_metadata_20190910.csv")

# Load the data.table that contains all of the somamer ratios as well as the identifier
soma_dt <- readRDS("./data/ms_ratios_20190910.rds")

# In the following code, we take the metadata, filter down to the training cohort, select only the
# Identifier and the trait of interest (msdss) and convert to a data.table to merge with the somamer ratios. 
# If we wanted to build a model predicting a different outcome with a different subset of data, those steps can
# be included here
data <- data %>% 
  filter(model_cohort=="training") %>% 
  select(sampleid,msdss) %>% 
  as.data.table()

# Next lines merge the retained metadata with the ratios and then removes the identifier.
# Note that the ranger R package uses all of the avaiable columns (other than the specified response) 
# as predictors. This is why the idenfifier is removed in the pipeline.
soma_dt <- merge(data,soma_dt,by="sampleid",sort=TRUE)
soma_dt <- soma_dt[,-1]

# Create a character vector of all somamer ratios (and individual markers) as well as outcomes of interest 
# for use downstream
all_rats <- names(soma_dt)

# Command line input using the batch R package. Tells R where to look for the ratios that we 
# want to include at current iteration. Then filters the data that will be used for modeling to 
# only include the outcome of interest and the ratios (somamers) at the current iteration
var_path <- "./scripts/msdss/iter0_ratios.txt"
parseCommandArgs()
if(var_path != ""){
	var_list <- readLines(var_path)
	soma_dt <- subset(soma_dt,select=which(all_rats %in% c("msdss",var_list)))
 }

# Selecting tuning parameter based on the size of the data
# if() statement was added to account for iterations with very few ratios
soma_dt <- as.data.frame(soma_dt)
nvars <- floor(3*sqrt(dim(soma_dt)[2]))
if(nvars >= dim(soma_dt)[2]-1){nvars <- floor(sqrt(dim(soma_dt)[2]))}

# DUMX is a required parameter for the Rswarm utility, but not needed here. Serves as a placeholder
hold <- DUMX

# Calculate the number of available CPUs for the forest, which is passed into the call to ranger
ncpus <- detectBatchCPUs()
mod <- ranger(data=soma_dt, # The full dataset containing predictors and response variables
              num.trees=40000, # number of trees in the forest
              mtry=nvars,      # number of variables to sample for each split
              importance='impurity',  # Type of variable importance measure to calculate
              num.threads=ncpus,      # Number of threads for parallelization
	            seed=DUMZ,              # setting a seed for reproducibility, DUMZ argument from Rswarm utility
	            dependent.variable.name="msdss")  # Telling ranger what the response variable is

saveRDS(mod,"DUMY1")  # Saving the model as an ouput, DUMY1 argument from Rswarm utility
