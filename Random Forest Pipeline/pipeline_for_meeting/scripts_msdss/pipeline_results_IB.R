library(readr, quietly = TRUE)
library(readxl, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(tidyr, quietly = TRUE)
library(lubridate, quietly = TRUE)
library(bindrcpp, quietly = TRUE)
library(stringr, quietly = TRUE)
library(ranger)
library(ggplot2)

install.packages("readr")
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("lubridate")
install.packages("bindrcapp")
install.packages("stringr")
install.packages("ranger")


source("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\Bibi\\pipeline_for_meeting\\pipeline_for_meeting\\scripts\\pipeline_results_function_IB.R")

setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\Bibi\\pipeline_for_meeting\\pipeline_for_meeting\\data")

soma_untreated <- read.csv("ms_metadata_20190910.csv")

ratio_append <- function(data,ratio_list){
  splitlist <- strsplit(ratio_list,"[/]")
  v1 <- sapply(splitlist,first)
  v2 <- sapply(splitlist,last)
  p <- length(ratio_list)
  for(i in 1:p){
    data[,ratio_list[i]] <- data[[v1[i]]]/data[[v2[i]]]
  }
  return(data)
}

setwd("C:\\Users\\bishofij\\Desktop\\here\\msdss")

msdss_plot <- pipeline_results("msdss","MS-DSS Baseline Visit",iters=0:119)
msdss_plot <- pipeline_results("msdss","MS-DSS Baseline Visit",iters=0:119,
                       best_iter=103)
msdss_plot

# Using best MSE
# msdss_new_iter <- 92

# Using eye ball
msdss_iter <- 103


msdss_var <- readLines(paste('C:\\Users\\bishofij\\Desktop\\here\\msdss\\iter',msdss_iter,"_ratios.txt",sep=""))
soma_untreated <- soma_untreated %>% ratio_append(msdss_var)

soma_train <- soma_untreated %>% filter(model_cohort=="training")
soma_test <- soma_untreated %>% filter(model_cohort=="validation")

nvars <- floor(3*sqrt(length(msdss_var)))
msdss_mod <- ranger(data=as.data.frame(soma_train[,c("msdss",msdss_var)]), num.trees=40000,mtry=nvars,importance='impurity',
              seed=223,dependent.variable.name="msdss",write.forest = TRUE)

######################################################################
##################   Add predictions to Dataset ######################
######################################################################
soma_train$soma_msdss <- msdss_mod$predictions
soma_test$soma_msdss <- predict(msdss_mod,soma_test)$predictions
soma_untreated <- bind_rows(soma_train,soma_test)

############################################################################
############### Ratios and Variable Importance in Models ###################
############################################################################
msdss_imp_dat <- tibble(ratio = names(msdss_mod$variable.importance),
                        importance = as.numeric(msdss_mod$variable.importance))

################################################################################
######################### Plot of observed versus predicted ####################
################################################################################
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

