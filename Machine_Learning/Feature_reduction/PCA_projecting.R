
# Project PCA
library("dplyr")

# Pull in data
setwd("~/data/horvath_data2")
data <- readRDS("age_balanced_qnorm_450k_GSE42861_bin10_age25.rds")
datMethTrain <- training[,-c(1:3)]

#Impute missing values if needed. You can use a different imputation method of your choice but we have not found 
#     this makes a significant difference.
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMethTrain <- apply(datMethTrain,2,meanimpute)


# Get test data
# Validation data
setwd("~/data/horvath_data2")
validation 

# Filter and prep Validation
validation <- filter(validation, study != "GSE42861")
validation <- filter(validation, between(age, 25, 70))
datMethTest <- validation[,-c(1:3)]
datMethTest <- apply(datMethTest,2,meanimpute)



#Perform PCA and projections. Remove last PC.
PCA = prcomp(datMethTrain,scale.=F)
TrainPCData = PCA$x[,1:(dim(PCA$x)[2]-1)]
TestPCData = predict(PCA,datMethTest)[,1:(dim(PCA$x)[2]-1)]

TrainPCData <- as.data.frame(TrainPCData)
TrainPCData$age <- training$age
TrainPCData$sample <- training$sample
training <- TrainPCData

TestPCData <- as.data.frame(TestPCData)
TestPCData$age <- validation$age
TestPCData$sample <- validation$sample
validation <- TestPCData


# PCA of cell lines and colored by indications
input <- data[,-c(1:2)]
input[is.na(input)] = 0
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
#title <- paste0("Test ",as.numeric(ncol(training)))

# Make PCA plot and print
library(ggplot2)
library(ggfortify)
plot <- autoplot(PCA,
                 data = training,
                 colour = 'age',
                 label = FALSE,
                 label.size = 4,
                 shape = 19,
                 #frame = TRUE,
                 main = "PCA")
print(plot)
plot(training$age, training$PC1)
