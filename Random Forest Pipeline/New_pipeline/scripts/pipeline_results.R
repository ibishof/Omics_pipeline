####################################################################################
# This script can be used to find the iteration of the pipeline that we want to use
# as the final model, recreate the model, add the model predictions to the metadata, 
# and create variable importance measures for the final model
####################################################################################
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(bindrcpp)
library(stringr)
library(ranger)
library(ggplot2)

# Load function that does the heavy lifting
source("./scripts/pipeline_results_function.R")

# Function that takes a list of ratios, appends them to the end of a dataset, 
# and returns the dataset
ratio_append <- function(dat,ratio_list){
  splitlist <- strsplit(ratio_list,"[/]")
  v1 <- sapply(splitlist,first)
  v2 <- sapply(splitlist,last)
  p <- length(ratio_list)
  for(i in 1:p){
    dat[,ratio_list[i]] <- dat[[v1[i]]]/dat[[v2[i]]]
  }
  return(dat)
}

# Load the metadata
soma_untreated <- read_csv("./data/ms_metadata_20190910.csv")

# Example of how to examine results from one pipeline, prints iteration
# with best training performance, and puts line at that point on plot
msdss_plot <- pipeline_results("msdss","MS-DSS Baseline Visit",iters=0:119)
msdss_plot

# Specifing best_iter places line at that location
msdss_plot <- pipeline_results("msdss","MS-DSS Baseline Visit",iters=0:119,
                       best_iter=109)
msdss_plot

# Pick the iteration that will be used for the final model developed
msdss_iter <- 103

# Loads in the set of ratios/somamers that were used at the specified iteration,
# then append them to the metadata
msdss_var <- readLines(paste('./scripts/msdss/iter',msdss_iter,"_ratios.txt",sep=""))
soma_untreated <- soma_untreated %>% ratio_append(msdss_var)

# Filter the data into training and validation cohorts
soma_train <- soma_untreated %>% 
  filter(model_cohort=="training")
soma_test <- soma_untreated %>% 
  filter(model_cohort=="validation")

# Build the model using ranger, using only the training data
nvars <- floor(3*sqrt(length(msdss_var)))
msdss_mod <- ranger(data=as.data.frame(soma_train[,c("msdss",msdss_var)]), 
                    num.trees=40000,mtry=nvars,importance='impurity',
              seed=223,dependent.variable.name="msdss",write.forest = TRUE)

# Save the model object for later use
saveRDS(msdss_mod,"./models/msdss_mod.rds")

# Can be loaded with the following:
# msdss_mod <- readRDS("./models/msdss_mod.rds")

######################################################################
##################   Add predictions to Dataset ######################
######################################################################

# For the training data, using the out-of-bag predictions from the model object
soma_train$soma_msdss <- msdss_mod$predictions

# For the validation data, create the predictions using the predict function as so
soma_test$soma_msdss <- predict(msdss_mod,soma_test)$predictions

# Merge the two back together
soma_untreated <- bind_rows(soma_train,soma_test)

############################################################################
############### Ratios and Variable Importance in Models ###################
############################################################################

# Create a dataset that contains the ratios that were used in the final model
# along with the variable importance measures for said ratios
msdss_imp_dat <- tibble(ratio = names(msdss_mod$variable.importance),
                        importance = as.numeric(msdss_mod$variable.importance))
write.csv(msdss_imp_dat, "msdss_imp_dat.csv")

################################################################################
######################### Plot of observed versus predicted ####################
################################################################################

# Now that you have all of the pieces, any comparisons can be done. In the following,
# I examine the correlation between the model predicted values versus the observed 
# values in the validation cohort

my_cor <- function(x,y,...){
  obj <- cor.test(x,y,...)
  val <- obj$estimate
  pval <- obj$p.value
  out <- paste(round(val,3)," (",format.pval(pval,eps=0.001,digits=4),")",sep="")
}
my_scatterplot <- function(x,y,main="",fix=TRUE,xlab=NULL,ylab=NULL,cor_method="spearman",reg=FALSE,
                           xloc = 3.5, yloc = .8, ...){
  if(fix==TRUE){
    if(is.null(xlab)){xlab <- "Predicted"}
    if(is.null(ylab)){ylab <- "Observed"}
    lim <- range(x,y)
    plot(x,y,xlim=lim,ylim=lim,xlab=xlab,ylab=ylab,
         main=main,cex.lab=1.5,cex.main=1.5,cex.axis=1.3,...)
    text(xloc,yloc,paste("Spearman =",my_cor(x,y,method=cor_method)),cex=1.5)
    if(reg==TRUE){
      abline(lm(y~x))
    }else{abline(0,1)}
  }  
  if(fix==FALSE){
    plot(x,y,xlab=xlab,ylab=ylab,
         main=main,cex.lab=1.5,cex.main=1.5,...)
    legend("topleft",paste("Spearman =",my_cor(x,y,method=cor_method)),cex=2,bty="n")
    text(xloc,yloc,paste("Spearman =",my_cor(x,y,method=cor_method)),cex=1.5)
    if(reg==TRUE){
      abline(lm(y~x))
    }else{abline(0,1)}}
}
with(soma_test,my_scatterplot(soma_msdss,msdss,main="MS-DSS Baseline Visit"))
