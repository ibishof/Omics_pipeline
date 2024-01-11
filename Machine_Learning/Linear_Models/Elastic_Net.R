


# Elastic net pipeline


library(caret)
library(glmnet)
library(dplyr)

# Optional select features to keep
setwd("~/data/horvath_data2/blood/modeling/effect_of_study/GSE42861A_qnorm_no_colinear")
top_features <- read.csv("271_variable_importance_calc.csv")
top_features <- top_features[order(-top_features$Rmean), , drop = FALSE]
features_to_keep <- top_features$Features

# training data
setwd("~/data/horvath_data2")
training <- readRDS("age_balanced_qnorm_450k_GSE42861_bin10_age25.rds")
traning_samples <- training$sample

data <- readRDS("age_balanced_qnorm_450k_GSE42861_bin10_age25.rds")
data <- filter(data, sample %in% traning_samples)
data <- data %>% select(where(is.numeric))
data <- select(data, age, top_cpg)
# data$rf_age <- training_predictions$predictions
x <- as.matrix(select(data, -age))

# Response Variable
y <- data$age


# Validation data
setwd("~/data/horvath_data2")
validation <- readRDS("qnorm_blood_450k.rds")

# Add metadata if needed
metadata <- readRDS("qnorm_blood_450k.rds")
metadata <- select(metadata, sample, study)
validation <- merge(metadata, validation, by = "sample")

# Filter and prep Validation
validation <- filter(validation, sample %in% validation_samples)
# validation$rf_age <- predictions$predictions

validation <- filter(validation, !(sample %in% traning_samples))
validation <- filter(validation, study != "GSE42861")
validation <- filter(validation, between(age, 25, 70))
# validation <- select(validation, age, features_to_keep)
# validation <- select(validation, age, top_cpg)
validation <- validation %>% select(where(is.numeric))
newx <- as.matrix(select(validation, -age))

# Cross-validation
cvfit <- cv.glmnet(x, y)

# Plots the cross-validation curve (red dotted line) along with upper and lower standard deviation curves along the λ sequence
plot(cvfit)

# Lamda at min mean sq error
cvfit$lambda.min
# Lambda that gives the most regularized model such that the cross-validated error is within one standard error of the minimum
cvfit$lambda.1se

# Coef at min sq error
coef(cvfit, s = "lambda.min") # Can set S to anything to explore other lamda values
# Make nice table
coef_lambda.min <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")))

# Make prediction with new model
predictions <- predict(cvfit, newx = x, s = "lambda.min")

# Plot predictions 
corr <- cor(data$age, predictions)
corr
plot(data$age, predictions, main = paste0("Training Corr = ", round(corr, 3)))
abline(coef = c(0,1), col = "red")

# Plot Validation prediction
# Make prediction with new model
predictions <- predict(cvfit, newx = newx, s = "lambda.min")

# Plot predictions 
corr <- cor(validation$age, predictions, use = "complete.obs")
corr
plot(validation$age, predictions, main = paste0("Validation Corr = ", round(corr, 3)))
abline(coef = c(0,1), col = "red")


median_error <- median(abs(validation$age - predictions), na.rm = TRUE)
plot(validation$age, (validation$age - predictions), 
     main = paste0("Error Validation, Median Error = ", round(median_error, 3)),
     xlab = "Age",
     ylab = "Error",
     pch = 1)
abline(h = 0, col = "red")

# Create table of non-zero coefficients
coef_lambda.min.nozero <- filter(coef_lambda.min, s1 != 0)
data_top <- select(data, age, rownames(coef_lambda.min.nozero)[-1])
rownames(data_top) <- data$sample


# 100x
coef_lambda.min.summary <- list()
median_error <- list()
for (i in 1:10){
  cvfit <- cv.glmnet(x, y, alpha = 0.5, lower.limits = -5, upper.limits = 5)
  # Make nice table
  coef_lambda.min.summary[[i]] <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
  
  # Calculate median error
  predictions <- predict(cvfit, newx = newx, s = "lambda.min")
  median_error[i] <- median(abs(validation$age - predictions), na.rm = TRUE)
  }
  
coef_lambda.min.table <- as.data.frame(coef_lambda.min.summary)
coef_lambda.min.table$mean <- rowMeans(coef_lambda.min.table)
coef_lambda.min.table <- filter(coef_lambda.min.table, mean != 0)

mean(unlist(median_error))

# Prepare data for ratio building
data <- readRDS("qnorm_blood_450k.rds")
sample <- data$sample
age <- data$age
data <- select(data, rownames(coef_lambda.min.table)[-1])
top_cpg <- rownames(coef_lambda.min.table)[-1]
