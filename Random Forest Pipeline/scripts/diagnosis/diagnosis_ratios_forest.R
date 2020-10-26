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

data <- read_csv("./data/ms_metadata_20190910.csv")
soma_dt <- readRDS("./data/ms_ratios_20190910.rds")

names(data)[names(data) == "Diagnosis_Main"] <- "diagnosis"
names(data)[names(data) == "SMF_ID_or_STS2_Identifier"] <- "sampleid"
names(soma_dt)[names(soma_dt) == "SMF_ID_or_STS2_Identifier"] <- "sampleid"
data$diagnosis[data$diagnosis != "MS"] <- "Non-MS"

data <- data %>% 
  filter(model_cohort=="training") %>% 
  select(sampleid,diagnosis) %>% 
  as.data.table()

soma_dt <- merge(data,soma_dt,by="sampleid",sort=TRUE)
soma_dt <- soma_dt[,-1]

all_rats <- names(soma_dt)

# INPUTS GO HERE
remove_outliers <- "FALSE"
var_path <- "./scripts/diagnosis/iter0_ratios.txt"
parseCommandArgs()
if(var_path != ""){
  var_list <- readLines(var_path)
  soma_dt <- subset(soma_dt,select=which(all_rats %in% c("diagnosis",var_list)))
}

# Selecting tuning parameter
soma_dt <- as.data.frame(soma_dt)
nvars <- floor(3*sqrt(dim(soma_dt)[2]))
if(nvars >= dim(soma_dt)[2]-1){nvars <- floor(sqrt(dim(soma_dt)[2]))}
hold <- DUMX

ncpus <- detectBatchCPUs()
soma_dt$diagnosis <- as.factor(soma_dt$diagnosis)
mod <- ranger(data=soma_dt, num.trees=min(40000,dim(soma_dt)[2]),mtry=nvars,importance='impurity',num.threads=ncpus,
              seed=DUMZ,dependent.variable.name="diagnosis",classification = TRUE,keep.inbag = TRUE)

sonar.task = makeClassifTask(data = soma_dt, target = "diagnosis")
results = OOBCurve(mod, task = sonar.task, data = soma_dt)
mod$oobauc <- results$auc[length(results$auc)]

saveRDS(mod,"DUMY1")
