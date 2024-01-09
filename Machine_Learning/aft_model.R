
# Build Accelerated Failure Time model (AFT)

# Load the packages
library(survival)
library(survminer)
library(dplyr)

# Fit the Accelerated Failure Time model
aft_weibull <- survreg(surv_object_train ~ ., 
                     data = select(training_set, -follow_up_time, -`Date G30 first reported (alzheimer's disease)`, -eid, -event),
                     dist = 'weibull')

aft_lognormal <- survreg(surv_object_train ~ ., 
                         data = select(training_set, -follow_up_time, -`Date G30 first reported (alzheimer's disease)`, -eid, -event), 
                         dist = 'lognormal')

aft_exponential <- survreg(surv_object_train ~ ., 
                           data = select(training_set, -follow_up_time, -`Date G30 first reported (alzheimer's disease)`, -eid, -event), 
                           dist = 'exponential')

aft_logistic <- survreg(surv_object_train ~ ., 
                           data = select(training_set, -follow_up_time, -`Date G30 first reported (alzheimer's disease)`, -eid, -event), 
                           dist = 'logistic')

aft_gaussian <- survreg(surv_object_train ~ ., 
                        data = select(training_set, -follow_up_time, -`Date G30 first reported (alzheimer's disease)`, -eid, -event), 
                        dist = 'gaussian')


# Check goodness of fit
AIC(aft_model,
    aft_weibull,
    aft_exponential,
    aft_lognormal,
    aft_logistic,
    aft_gaussian)
saveRDS(aft_weibull, "aft_weibull.rds")
saveRDS(aft_exponential, "aft_exponential.rds")
saveRDS(aft_exponential, "aft_exponential.rds")
saveRDS(aft_lognormal, "aft_lognormal.rds")
saveRDS(aft_gaussian, "aft_gaussian.rds")

summary(aft_model)

# Extract residuals
dev_res <- residuals(aft_model, type = "deviance")
# Get the theoretical Weibull quantiles
weibull_quantiles <- qweibull(ppoints(length(dev_res)), shape = 1, scale = 1)

# Create Q-Q plot
qqplot(weibull_quantiles, sort(dev_res))
abline(0, 1) # Adds 45-degree line



# Build simple model using just age and sex
simpler_model <- survreg(surv_object ~ age_at_blood_draw + sex, data = data)


# Perform likelihood ratio test
lr_test <- anova(simpler_model, aft_model)
print(lr_test)

# Comapre to cox model
logLik(cox_model)
logLik(aft_model)

AIC(aft_model)
AIC(cox_model)

coeff_table <-  as.data.frame(aft_model$coefficients[-1])
sd1 <- apply(select(training_set, -follow_up_time, -`Date G30 first reported (alzheimer's disease)`, -eid, -event), 2, sd)
coeff_table$importance_sd <- sd1*coeff_table$`aft_model$coefficients`


# Predictions


# Making predictions for specific observations
predictions <- predict(aft_model, newdata = validation_set)
quantile(predictions, probs = seq(0,1,.01))
hist(predictions)
min(predictions)
max(predictions)


cox_surv_fit <- survfit(cox_model, newdata = validation_set)
cox_probabilities <- summary(cox_surv_fit, times = c(1:10))$surv


# Get c-index
cox_pred <- predict(cox_model, newdata = validation_set, type = "risk")
aft_model <- aft_weibull
aft_pred <- exp(predict(aft_model, newdata = validation_set))
# Replace inf with max non-inf value
max_val <- max(aft_pred[aft_pred != Inf], na.rm = TRUE)
aft_pred[is.infinite(aft_pred)] <- max_val

# Assuming `time` and `status` are the columns in your validation_set
library(Hmisc)
c_index_cox <- rcorr.cens(cox_pred, Surv(validation_set$`Date G30 first reported (alzheimer's disease)`,
                                         validation_set$event))

c_index_aft <- rcorr.cens(aft_pred, Surv(validation_set$`Date G30 first reported (alzheimer's disease)`,
                                         validation_set$event))
