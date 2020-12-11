# Isolate QC and QA samples
library(data.table)

setwd("C:\\Users\\kqjf682\\Omics_pipelines\\msCleanr\\full_batch")
data1 <- read.csv("Height_0__lowessIS_20201241540.csv", row.names = "Alignment.ID", stringsAsFactors = FALSE)
datat1 <- as.data.frame(t(data1))
QC_height <- (datat1[rownames(datat1) %like% "QC", ])
QA_height <- (datat1[rownames(datat1) %like% "QA", ])

setwd("C:\\Users\\kqjf682\\Omics_pipelines\\msCleanr\\full_batch")
data2 <- read.csv("Normalized_0_lowessIS_20201241540.csv", row.names = "Alignment.ID", stringsAsFactors = FALSE)
datat2 <- as.data.frame(t(data2))
QC_Normalized <- (datat2[rownames(datat2) %like% "QC", ])
QA_Normalized <- (datat2[rownames(datat2) %like% "QA", ])

setwd("C:\\Users\\kqjf682\\Omics_pipelines\\msCleanr\\full_batch\\final_data")
MScleanR_QC <- read.csv("annotated_MS_peaks-cleaned_QC.csv")

# QA_height$batch <- MScleanR_QC$Batch
# QA_Normalized$batch <- MScleanR_QC$Batch


setwd("C:\\Users\\kqjf682\\Omics_pipelines\\msCleanr\\full_batch\\final_data")
metadata <- read.csv("metadata.csv")

QA_QC_height <- rbind(QC_height, QA_height)
QA_QC_normalized <- rbind(QC_Normalized, QA_Normalized )

QA_QC_height$batch <- metadata$Batch
QA_QC_normalized$batch <- metadata$Batch

##################################################################
library(BBmisc)
QA_QC_heightz <- normalize(QA_QC_height, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
QA_QC_heightz$batch <- metadata$Batch

for (col in 1:ncol(QA_QC_heightz)) {
  if (exists("htable")){
    if (ncol(QA_QC_heightz) != col) {
    itable <- cbind(metadata$Batch, QA_QC_heightz[,col])
    htable <- rbind(htable, itable)}
  }
  if (isFALSE(exists("htable"))){
    htable <- QA_QC_heightz[,col]
    htable <- cbind(metadata$Batch, htable)
  }
 
}
htable <- as.data.frame(htable)

htable <- htable %>% 
  rename(
    batch = V1,
    zscores = htable
  )


QA_QC_normalizedz <- normalize(QA_QC_normalized, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
QA_QC_normalizedz$batch <- metadata$Batch

for (col in 1:ncol(QA_QC_normalizedz)) {
  if (exists("ntable")){
    if (ncol(QA_QC_normalizedz) != col) {
      jtable <- cbind(metadata$Batch, QA_QC_normalizedz[,col])
      ntable <- rbind(ntable, jtable)}
  }
  if (isFALSE(exists("ntable"))){
    ntable <- QA_QC_normalizedz[,col]
    ntable <- cbind(metadata$Batch, ntable)
  }
  
}
ntable <- as.data.frame(ntable)

ntable <- ntable %>% 
  rename(
    batch = V1,
    zscores = ntable
  )
##################################################################

library(ggplot2)
# Basic box plot
htable$batch <- as.factor(htable$batch)
p <- ggplot(htable, aes(x=batch, y=zscores, fill = batch)) + 
  ggtitle("Height Range") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot()
p

ntable$batch <- as.factor(ntable$batch)
p <- ggplot(ntable, aes(x=batch, y=zscores, fill = batch)) + 
  ggtitle("Normalized Range") + theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot()
p



