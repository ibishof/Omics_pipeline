

# Remove highly correlated features
tmp <- cor(training)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
# Above two commands can be replaced with 
# tmp[!lower.tri(tmp)] <- 0
#
 
   data.new <- data[,!apply(tmp,2,function(x) any(x > 0.30))]
   head(data.new)
  
   
   library(caret)
   df2 = cor(training)
   hc = findCorrelation(df2, cutoff=0.10) # putt any value as a "cutoff" 
   hc = sort(hc)
   reduced_Data = df1[,-c(hc)]
   print (reduced_Data)