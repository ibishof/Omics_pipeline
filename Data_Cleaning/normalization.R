
library(matrixStats)
library(BBmisc)
library(dplyr)
library(tidyverse)

# Read in data
setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\module_reduced\\module_reduced_CKD_Diabetes\\svm")
train <- read.csv("training.csv")
train <- train[,-1]
validation <- read.csv("test.csv")
validation <-validation[,-1]

#####Various normalization methods
#
# Z-score # 
# must be numeric only
trainN <- train[2:ncol(train)]
train_mean <- colMeans(trainN)
train_mean <- sweep(trainN, 2, train_mean, FUN = '-')
trainsd <- t(as.matrix(colSds(as.matrix(trainN))))
trainz <- sweep(train_mean, 2, trainsd, FUN = '/')
train <- cbind(train$Diagnosis, data.frame(trainz))

validationN <- validation[2:ncol(validation)]
validation_mean <- colMeans(validationN)
validation_mean <- sweep(validationN, 2, validation_mean, FUN = '-')
validationsd <- t(as.matrix(colSds(as.matrix(validationN))))
validationz <- sweep(validation_mean, 2, validationsd, FUN = '/')
validation <- cbind(validation$Diagnosis, data.frame(validationz))

train1 <- train %>% 
  rename(
    Diagnosis = `train$Diagnosis`
  )

validation1 <- validation %>% 
  rename(
    Diagnosis = `validation$Diagnosis`
  )

# Same Z-score just less code
# Need library(BBmisc)
library(BBmisc)
train <- normalize(train, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")
validation <- normalize(validation, method = "standardize", range = c(0, 1), margin = 1L, on.constant = "quiet")

# Range from 0 to 1
# Normalizing method. Available are: "center": Subtract mean.
# "scale": Divide by standard deviation. 
# "standardize": Center and scale. "range": 
# Scale to a given range.
#
train <- normalize(train, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
validation <- normalize(validation, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")


# Normalize by dividing by the max value of each column
#
trainN <- train[2:ncol(train)]
train_max <- apply(trainN,2,max)
train_norm <- sweep(trainN, 2, train_max, FUN="/")

validationN <- validation[2:ncol(validation)]
validation_max <- apply(validationN,2,max)
validation_norm <- sweep(validationN, 2, validation_max, FUN="/")
