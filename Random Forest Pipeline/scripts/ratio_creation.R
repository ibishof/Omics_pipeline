##################################################################
# Script that can be used to create a data.table containing 
# all possible somamer ratios, in addition to individual somamers
##################################################################

# Libraries that are used
library(readr)
library(dplyr)
library(tidyr)
library(data.table)
library(foreach)

# Function that takes a character ratio name and splits into numerator and denompnator
split_rat <- function(rat){
  split_list <- unlist(strsplit(rat,"[/]"))
  v1 <- first(split_list)
  v2 <- last(split_list)
  return(c(v1,v2))
}

# Creates a ratio of somamers (rat) given that they somamers are in the dataset (dat)
create_rat <- function(rat,dat){
  ind_mark <- split_rat(rat)
  out <- dat[[ind_mark[1]]]/dat[[ind_mark[2]]]
  return(out)
}

# Load the metadata
soma_untreated <- read_csv("/Users/wangq13/Desktop/Bibi/adjusted_metadata_claire.csv")
highsignal_markers <- read_csv("/Users/wangq13/Desktop/Bibi/top_1000_somamers.csv", col_types = cols())
highsignal_markers <- as.data.frame(highsignal_markers)

# Create a character vector of all the Somamers in the dataset
# NOTE: THE SL CODES ARE NOT UNIQUE IN THE 5K DATA, SO THIS LINE WILL NEED TO BE 
# ALTERED WHEN TRANSITIONING TO THE 5K DATA.

markers <- highsignal_markers$x

# Reduce the data to only include the unique identifier (sampleid in this case) and 
# all the individual Somamers.
soma_untreated <- soma_untreated %>% 
  select(SMF_ID_or_STS2_Identifier,which(names(soma_untreated) %in% markers))

# Create a character vector of all possible ratio combinations
ratios <- character(length=choose(length(markers),2))
iter <- 1
for(i in 1:(length(markers)-1)){
  for(j in (i+1):length(markers)){
    ratios[iter] <- paste(markers[i],"/",markers[j],sep="")
    iter <- iter + 1
  }
}

# Create a data.table of all possible ratio combinations using the foreach R package
# Allows for parallel implementations using the %dopar% function and proper backendd
soma_dt <- foreach(i=1:length(ratios),.final=as.data.table,.packages=c('dplyr')) %do% {
  create_rat(ratios[i],dat=soma_untreated) 
  } 
names(soma_dt) <- ratios

# Assign unique identifiers
soma_dt[,SMF_ID_or_STS2_Identifier := soma_untreated$SMF_ID_or_STS2_Identifier]
setcolorder(soma_dt,c("SMF_ID_or_STS2_Identifier"))

# Add individual somamers to data.table, merge with ratios, and save
markers_dt <- soma_untreated %>% 
  as.data.table()
soma_dt <- merge(soma_dt,markers_dt,by="SMF_ID_or_STS2_Identifier")

temp_name <- gsub("/", "_", names(soma_dt))
names(soma_dt) <- temp_name

saveRDS(soma_dt,"/Users/wangq13/Desktop/Bibi/ms_ratios_20190910.rds")

# Creates a character vector of ratios and individual somamers, can convert this to "iter0_ratios.txt" for all of the
# pipelines that will be run.
ratios <- gsub("/", "_", ratios)
writeLines(c(ratios,markers),"/Users/wangq13/Desktop/Bibi/iter0_ratios.txt")
