library(e1071)
library(matrixStats)
library(BBmisc)
library(psycho)
library(dplyr)
library(tidyverse)
help(package="e1071")


setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\module_reduced\\module_reduced_CKD_Diabetes\\svm")
train <- read.csv("training.csv")
train <- train[,-1]
validation <- read.csv("test.csv")
validation <-validation[,-1]

# Normalize
#
train <- normalize(train, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")
validation <- normalize(validation, method = "range", range = c(0, 1), margin = 1L, on.constant = "quiet")

# Scale data by calculating z-score
# 
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

 train <- train %>% 
  rename(
    Diagnosis = `train$Diagnosis`
      )
 
 validation <- validation %>% 
   rename(
     Diagnosis = `validation$Diagnosis`
   )
 


# Build the actual model
# The smaller gamma is, the more the hyperplane is going to look like a straight line. If gamma is too great, the hyperplane will be more curvy 
# cost = the penaility of constraints violation
# eps-regression is default for regression models, nu-regression can be selected also
# (default) C-classification, nu-classification, one-classification (for novelty detection)
# Kernal defualt is radial; linear, polynomial, and sigmoid are also options
model <- svm(Diagnosis ~ ., data = train, gamma = 1e-05, cost = .01, type = "C-classification", kernel = "linear")
summary(model)

# Test model on training regression
pred <- predict(model, train)
cor(pred,soma_dt$GFR_0, method = "pearson")
table(pred, y)


# Test model for training classification
pred <- predict(model, train)
train_pred <- as.factor(pred)
tyyyyy <- ordered(train_pred)
roc_ojb_train <- roc(response =  train$Diagnosis, predictor =tyyyyy )
roc_ojb_train$auc

# Test model on validation
#validation <- read.csv("test.csv")
pred <- predict(model, validation)
cor(pred,validation$GFR_0, method = "pearson")
table(pred, y)

#Get AUC
pred <- predict(model, validation)
zzzz <- as.factor(pred)
yyyyy <- ordered(zzzz)
roc_ojb_MS <- roc(response =  validation$Diagnosis, predictor =yyyyy )
roc_ojb_MS$auc

write.csv(pred, "svm_MS.csv")
write.csv(answers, "answers-svm_ms.csv")


# Tune svm parameters
tuned_parameters <- tune.svm(Diagnosis~., data = train, kernel = "linear", gamma = 10^(-5:-1), cost = 10^(-3:1))
summary <-summary(tuned_parameters )
summary

# COde to gain feature importance (Needs work)
library(rminer)
M <- fit(y~., data=soma_dt, model, kpar=list(sigma=0.10), C=2)
svm.imp <- Importance(M, data=soma_dt)



