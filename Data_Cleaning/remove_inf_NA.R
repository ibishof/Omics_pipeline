
#This scriptremoves columsn with 50% or more zeros. Replaces zero to NA. Log2 tranforms. Then converts the na back to zero.
setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\third")

soma_dt <-read.csv("random_forest_input.csv", row.names = 1)

soma_dt<- readRDS("Ms1_ratio.rds")
# Replace Inf with NA
for (j in 1:ncol(soma_dt)) set(soma_dt, which(is.infinite(soma_dt[[j]])), j, NA)

# Change NA to zero
soma_dt <- as.data.frame(soma_dt)
soma_dt <- soma_dt %>% mutate_if(is.numeric , replace_na, replace = 0)


