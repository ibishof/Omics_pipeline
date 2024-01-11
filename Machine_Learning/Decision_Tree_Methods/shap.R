



# Install the packages if you haven't already
if (!requireNamespace("iml", quietly = TRUE)) {
  install.packages("iml")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

setwd("~/scripts/scratch")
source("pablo_shap.R")

# Load the required packages
library(iml)
library(ggplot2)

data <- select(data, Acuity_0, X32_variable_importance_calc$Features)
x <- select(training, -Acuity_0)

pfun <- function(object, newdata) {
  # Get the probability predictions from the model
  prob_preds <- predict(object, data = newdata)$predictions
  
  # Convert factor predictions to numeric probabilities
  if (is.factor(prob_preds)) {
    prob_preds <- as.numeric(prob_preds == "sicker")  # Replace "sicker" with the name of the positive class
  }
  
  return(prob_preds)
}


# Compute fast (approximate) Shapley values using 10 Monte Carlo repetitions
system.time({  # estimate run time
  set.seed(5038)
  shap <- fastshap::explain(rating_mod, X = x, pred_wrapper = pfun, nsim = 100)
})


theme_set(theme_bw())

# Aggregate Shapley values
shap_imp <- data.frame(
  Variable = names(shap),
  Importance = apply(shap, MARGIN = 2, FUN = function(x) sum(abs(x)))
)



  
  # Plot Shap-based variable importance
  ggplot(shap_imp, aes(reorder(Variable, Importance), Importance)) +
  geom_col() +
  coord_flip() +
  xlab("") +
  ylab("mean(|Shapley value|)")

  
  # Make shap_result_bike with RF model
  # Order features by shap abs values
  ranks <- colMeans(abs(shap))
  shap <- shap[,order(-ranks)]
  shap_result_bike <- list()
  shap_result_bike$shap_score <- shap
  shap_result_bike$mean_shap_score <- colMeans(abs(shap))
  
  # Make bike_xr matrix
  bike_dmyr = dummyVars(" ~ .", data = x, fullRank=T)
  bike_xr = predict(bike_dmyr, newdata = x)
  
  
  
  ## Plot var importance based on SHAP
  var_importance(shap_result_bike, top_n=14)
  
  ## Prepare data for top N variables
  shap_long_bike = shap.prep(shap = shap_result_bike,
                             X_train = bike_xr , 
                             top_n = 14
  )
  
  
  
  
  ## Plot shap overall metrics
  plot.shap.summary(data_long = shap_long_bike)
  