
#This script removes columns with 50% or more zeros. Replaces zero to NA. Log2 tranforms. Then converts the na back to zero.
setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\WGCNA\\data_sets")

proteinGroups <-read.csv("72.csv", row.names = 1)

# Change zero to NA
soma_dt1[soma_dt1 == 0] <- NA

# Remove rows with 50% or more NA
delete.na <- function(DF, n=0) {
  DF[rowSums(is.na(DF)) <= n,]
}
cleandat <- delete.na(proteinGroups, ((ncol(proteinGroups))/2))

#Log 2 tranform and change NA back to zero
Log2 <-log2(cleandat)
Log2[is.na(Log2)] <- 0
table(is.na(Log2))

write.csv(Log2,"transformed.csv")


# Change inf to zero
validation[which(validation==Inf)] <- 0

# Replace NA with 0
validation[is.na(validation)] <- 0

