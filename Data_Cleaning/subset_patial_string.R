
setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\A003_1000")
data <- read.csv('Area_1_202010112021_unique_ms2_eGFR.csv')

# Subsetting table rows by partial string
library(data.table)
datat <- as.matrix(t(data))
standards <- (datat[rownames(datat) %like% "U13C", ])
write.csv(standards, "standards.csv")
data <- (standards)


library(stringi)
data1 <- datat[str_detect(datat, "U13"), ]  
