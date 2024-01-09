


# Elastic net pipeline


library(caret)
library(glmnet)
library(dplyr)
library(pROC)

# Prep data
setwd("~/data/MGH_Olink_COVID_Apr_27_2021")
data <- readRDS("npx_RF.rds")
data <- filter(data, COVID == 1)

# Create new feature
# data$change_in_severity <- data$Acuity_0 - data$Acuity_max
data <- select(data, -Acuity_max)
data <- data %>%
  mutate(binary_severity = ifelse(Acuity_0 %in% c(1,2), "sever", "mild"))
# data$change_in_severity <- as.factor(data$change_in_severity)
data$binary_severity <- as.factor(data$binary_severity)
data <- select(data, -Acuity_0)
# data <- select(data, -change_in_severity)

# Filter data
# Remove rows with too many NaN
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data[is.nan(data)] <- NA
data[is.na(data)] = 0

# Look for columns with charater variables
character_columns <- sapply(data, is.character)
character_column_names <- names(data)[character_columns]

# Get top features
top_features <- X19_variable_importance_calc$Features

# data$rf_age <- training_predictions$predictions
x <- as.matrix(select(data[,-1], top_features))

# Response Variable
y <- as.factor(data$binary_change)

# Encode the response variable as a binary outcome
y <- as.numeric(y) - 1


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
cvfit <- cv.glmnet(x, y, family="binomial", alpha=.01)

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

# Make predictions on the training data
predicted_probabilities <- predict(cvfit, x, type = "response")
predicted_classes <- ifelse(predicted_probabilities > 0.2, 1, 0)

# Calculate performance metrics
# Confusion matrix
confusion_matrix <- table(Predicted = predicted_classes, Actual = y)
print(confusion_matrix)

# Accuracy
accuracy <- mean(predicted_classes == y)
cat("Accuracy:", accuracy, "\n")

# Sensitivity and specificity
sensitivity <- confusion_matrix[2, 2] / (confusion_matrix[2, 2] + confusion_matrix[1, 2])
specificity <- confusion_matrix[1, 1] / (confusion_matrix[1, 1] + confusion_matrix[2, 1])
cat("Sensitivity:", sensitivity, "\n")
cat("Specificity:", specificity, "\n")

# Area under the curve (AUC)
roc_obj <- roc(as.numeric(y), as.numeric(predicted_probabilities))
auc_value <- roc_obj$auc
cat("AUC:", auc_value, "\n")

# Plot the ROC curve with the red line
roc_plot <- ggroc(roc_obj) + 
  labs(title = "Prediction of severity using proteins levels at day 0", x = "False Positive Rate", y = "True Positive Rate") +
  geom_abline(intercept = 1, slope = 1, linetype = "dotted", color = "black") +
  theme_bw()

roc_plot$layers[[1]]$aes_params$colour <- "red"
print(roc_plot)

# Boxplot
predicted_probabilities_df <- cbind(y, predicted_probabilities)

#
coef_lambda.min.table <- filter(coef_lambda.min, s1 != 0)


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
auc_value_list <- list()
for (i in 1:1000){
  indexes=sample(1:nrow(x),size = 1/5*nrow(x))
  newx=x[indexes,]
  y2=y[indexes]
  xtrain=x[-indexes,]
  ytrain=y[-indexes]
  
  cvfit <- cv.glmnet(xtrain, ytrain, family="binomial", alpha=.001)
  # Make nice table
  coef_lambda.min.summary[[i]] <- as.data.frame(as.matrix(coef(cvfit, s = "lambda.min")))
  
  # Calculate median error
  predicted_probabilities <- predict(cvfit, newx, type = "response")
  
  # Area under the curve (AUC)
  roc_obj <- roc(as.numeric(y2), as.numeric(predicted_probabilities))
  auc_value_list[i] <- as.numeric(sub("Area under the curve: ", "", roc_obj$auc))
}

coef_lambda.min.table <- as.data.frame(coef_lambda.min.summary)
coef_lambda.min.table$mean <- rowMeans(coef_lambda.min.table)
coef_lambda.min.table <- filter(coef_lambda.min.table, mean != 0)
coef_lambda.summary <- select(coef_lambda.min.table, mean)

mean(unlist(auc_value_list))

boxplot(unlist(auc_value_list), col = "orange", main = "AUC 1000 repeats")

# Prepare data for ratio building
data <- readRDS("qnorm_blood_450k.rds")
sample <- data$sample
age <- data$age
data <- select(data, rownames(coef_lambda.min.table)[-1])
top_cpg <- rownames(coef_lambda.min.table)[-1]






