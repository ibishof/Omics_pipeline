
library(xgboost)

# Bring training data
setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\module_reduced\\module_reduced_noC")
training <- read.csv("training.csv", check.names=FALSE)
training <- training[,-1]
training <- training[,order(colnames(training))]

#Bring in validation data
validation <- read.csv("test.csv", check.names=FALSE)
validation <- validation[,-1]
names <- unlist(colnames(training))
validation <-validation %>% 
  select(Diagnosis,names)
validation <- validation[,order(colnames(validation))]

# Convert the Diagnosis factor to an integer class starting at 0
# This is picky, but it's a requirement for XGBoost
state = training$Diagnosis
label = as.integer(state)-1
training$Diagnosis = NULL
training <- as.matrix(training)

state_v = validation$Diagnosis
label_v = as.integer(state_v)-1
validation$Diagnosis = NULL
validation <- as.matrix(validation)

# Transform the two data sets into xgb.Matrix
xgb.train = xgb.DMatrix(data=training,label=label)
xgb.test = xgb.DMatrix(data=validation,label=label_v)



# Parameters for grid search
num_class = length(levels(state))
params = list(
  booster="gblinear",
  eta=0.1,
  #max_depth=8,
  gamma=3,
  #subsample=0.75,
  colsample_bytree=1,
  objective="multi:softprob",
  eval_metric="merror",
  num_class=2,
  validate_parameters = TRUE
)

# Grid search for parameters max depth and child
max.depths = c(4:8)
child = c(1:4)

best_params = 1
best_score = 1

count = 1
for( depth in max.depths ){
  for( num in child){
    bst_grid=xgb.train(
      params=params,
      max.depth = depth,
      min_child_weight =num,
      #eta=num,
      data=xgb.train,
      nrounds=10000,
      nthreads=1,
      early_stopping_rounds=50,
      watchlist=list(val1=xgb.train),
      verbose=0
    )  
    
    if(count == 1){
      best_params = bst_grid$params
      best_score = bst_grid$best_score
      count = count + 1
    }
    else if( bst_grid$best_score < best_score){
      best_params = bst_grid$params
      best_score = bst_grid$best_score
    }
  }
}

best_params
best_score

# Parameters Second round
# Parameters for grid search
num_class = length(levels(state))
params = list(
  booster="gblinear",
  max_depth=best_params$max_depth,
  min_child_weight =best_params$min_child_weight,
  gamma=3,
  #subsample=0.75,
  colsample_bytree=1,
  objective="multi:softprob",
  eval_metric="merror",
  num_class=2,
  validate_parameters = TRUE
)


# Parameter search eta
etas = c(1, 0.1, 0.001, 0.0001, 0.00001)

best_params = 1
best_score = 1

count = 1
  for( num in etas){
    bst_grid=xgb.train(
      params=params,
      eta=num,
      data=xgb.train,
      nrounds=10000,
      nthreads=1,
      early_stopping_rounds=50,
      watchlist=list(val1=xgb.train),
      verbose=0
    )  
    
    if(count == 1){
      best_params = bst_grid$params
      best_score = bst_grid$best_score
      count = count + 1
    }
    else if( bst_grid$best_score < best_score){
      best_params = bst_grid$params
      best_score = bst_grid$best_score
    }
  }


best_params
best_score

# After tuning
# Use best parameters for model
num_class = length(levels(state))
params = list(
  booster= best_params$booster,
  eta= best_params$eta,
  max_depth=best_params$max_depth,
  gamma=best_params$gamma,
  #subsample=best_params$subsample,
  colsample_bytree=best_params$colsample_bytree,
  objective=best_params$objective,
  eval_metric=best_params$eval_metric,
  num_class=best_params$num_class,
  min_child_weight =3,
  validate_parameters = best_params$validate_parameters
)

# Train the XGBoost classifer
xgb.fit=xgb.train(
  params=params,
  data=xgb.train,
  nrounds=10000,
  nthreads=1,
  early_stopping_rounds=50,
  watchlist=list(val1=xgb.train),
  verbose=0
)

# Review the final model and results
xgb.fit
xgb.fit$best_score

# Predict outcomes with the test data
xgb.pred = predict(xgb.fit,validation,reshape=T)
xgb.pred = as.data.frame(xgb.pred)
colnames(xgb.pred) = levels(state)

#  identify the label (class) with the highest probability
xgb.pred$prediction = apply(xgb.pred,1,function(x) colnames(xgb.pred)[which.max(x)])
xgb.pred$label = levels(state)[label_v+1]

# Calculate the final accuracy
result = sum(xgb.pred$prediction==xgb.pred$label)/nrow(xgb.pred)
print(paste("Final Accuracy =",sprintf("%1.2f%%", 100*result)))


# Feature importance
var.importance <- xgb.importance(colnames(training), model = xgb.fit)
var.importance1 <- xgb.importance(colnames(training), model = xgb.fit, 
                                  data = xgb.train, label = label)

saveRDS(xgb.fit, "noC_boost.rds")
