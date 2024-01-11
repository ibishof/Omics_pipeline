

# This script loads training and validation datasets to fit an elastic net survival model
# using the glmnet package. It predicts risk scores for both datasets,
# calculates the model's c-index, and visualizes the average event value by decile of predicted scores.
# Feature importance is computed, and key results are saved for future analysis.



# Things you will need
# Training data
# Validation data (I use 10%)
# Event a binary vector of 0 or 1 indicating if the event of interest (AD diagnosis) has happened
# Your data will need a "time in years" that is length of follow up for no even sample. Sample with the event it time from baseline to that event.

setwd("your/path/here")
training_set <- read.csv("training_set.csv")
validation_set <- read.csv("validation_set.csv")

library(dplyr)

# Subsample 10% of the training_set
training_set <- training_set %>% sample_n(floor(0.1 * nrow(training_set)))
validation_set <- training_set %>% sample_n(floor(0.1 * nrow(validation_set)))


# Elastic net survival model
# The glmnet is a well documented elastic net model
library(glmnet)
library(survival)
library(dplyr)
library(Hmisc)


# Prepare the data
# Remove ids, length of follow up, and event status
# x_data is training
# y_data is your survival object with the time to event and event columns
# x2 is the validation data set
# time_diff_in_years
# Load required library

# Prepare the data
x_data <- as.matrix(select(training_set, -time_diff_in_years, -eid, -event, -X))
y_data <- as.matrix(Surv(training_set$time_diff_in_years, training_set$event))

# Fit the Lasso model
# alpha of alpha = 1: Lasso penalty. 
# Forces coefficients being exactly zero, effectively selecting a simpler model with fewer predictors.
# alpha = 0: Ridge regression.
# Ridge regression adds a penalty equivalent to the square of the magnitude of coefficients. All coefficients are shrunk by the same factor (none are eliminated).
fit <- cv.glmnet(x_data, y_data, family = "cox", alpha = .1)

# Prepare the validation set feature matrix
x_val <- as.matrix(select(validation_set, -time_diff_in_years, -eid, -event, -X))

# Predict risk scores for validation set
predicted_scores <- predict(fit, newx = x_val, type = "link", s = fit$lambda.min)
validation_set$predicted_scores <- predicted_scores
v_scores <- select(validation_set, eid, predicted_scores)

# Predict risk scores for training set
predicted_scores2 <- predict(fit, newx = x_data, type = "link", s = fit$lambda.min)
training_set$predicted_scores <- predicted_scores2
t_scores <- select(training_set, eid, predicted_scores)

# Create a Surv object for the validation set
y_val <- Surv(validation_set$time_diff_in_years, validation_set$event)


# Calculate the c-index
# This is similiar to AUC in that 1 is a perfect score and zero is the lowiest
c_index_result <- rcorr.cens(-predicted_scores, y_val)  # Notice the negative sign, lower risk score should be better
c_index_result
options(scipen=999)



# Visualizations
# Plot
plot(fit)


# Risk vs Age
plot(validation_set$age_at_blood_draw, validation_set$predicted_scores, col = rgb(0,0,1,alpha=0.1), pch = 19,
     xlab = "Age at Blood Drawn", ylab = "Hazard")
cor(validation_set$age_at_blood_draw, validation_set$predicted_scores)


# Plot the data
library(ggplot2)
# Sort validation_set by predicted_scores in descending order
validation_set <- validation_set[order(-validation_set$predicted_scores),]

# Add a column to indicate which decile each row belongs to
validation_set$decile <- cut(validation_set$predicted_scores, 
                              breaks = quantile(validation_set$predicted_scores, 
                                                 probs = 0:100/100),
                              labels = FALSE,
                              include.lowest = TRUE)

# Calculate the mean of 'event' for each decile
average_by_decile <- aggregate(event ~ decile, data = validation_set, FUN = mean)

# Plot the data as an XY plot
ggplot(average_by_decile, aes(x = as.factor(decile), y = event)) +
  geom_point() + # use geom_point for an XY plot
  geom_line(aes(group = 1)) + # connect the points with a line
  labs(x = "Decile of Predicted Scores",
       y = "Average Event Value",
       title = "Average Event Value by Decile of Predicted Scores")


basic_model <- average_by_decile
protein_model <- average_by_decile


# Now plot simple model vs fancy model
ggplot(basic_model, aes(x = as.factor(decile), y = event*100)) +
  geom_point(aes(color = "Age, Sex, PackYrs")) +  # First set of points
  geom_line(aes(group = 1, color = "Age, Sex, PackYrs")) +  # Connect the first set of points with a line
  geom_point(data = protein_model, aes(x = as.factor(decile), y = event*100, color = "Proteins")) +  # Second set of points
  geom_line(data = protein_model, aes(x = as.factor(decile), y = event*100, group = 1, color = "Proteins")) +  # Connect the second set of points with a line
  geom_hline(yintercept = (sum(validation_set$event)/nrow(validation_set))*100, linetype="solid", color = "black") +  # Add horizontal line
  labs(x = "Decile of Predicted Scores",
       y = "Average Event Value",
       title = "Average Event Value by Decile of Predicted Scores") +
  scale_x_discrete(breaks = seq(0, 100, by = 2)) +  # Set the tick marks on x-axis every 2
  scale_color_manual(values = c("Age, Sex, PackYrs" = "red", "Proteins" = "blue"))  # Set colors for the two sets of points



# Put feature importance into a single table
# Coef at min sq error
# Make nice table
coef_lambda.min <- as.data.frame(as.matrix(coef(fit, s = "lambda.min")))

coef_lambda.min$importance <- apply(x_data, 2, sd) * coef_lambda.min$`1`

saveRDS(fit, "enet.rds")
write.csv(coef_lambda.min, "coef_table.csv")
saveRDS(validation_set, "validation_set.rds")
saveRDS(training_set, "training_set.rds")
saveRDS(v_scores, "v_scores.rds")
saveRDS(t_scores, "t_scores.rds")
