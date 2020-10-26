# script to create correlation matrix plot
library(corrplot)
library(ggplot2)

pc <- clean_numbers

#define the columns that contain your abundance data. Change the number after the ":" to subset your data
com = pc#[,2:20]

# Now create a correlation matrix with your community composition data using the command 'cor':
cc = cor(com, method = "spearman")
###if you have missing values in your data add the 'use' parameter. Check '?cor' for the cor command help page
# I usually use Spearman correlation because I'm not overly concerned that my relationships fit a linear model,
#and Spearman captures all types of positive or negative relationships (i.e. exponential, logarithmic)


# Plot correlation
corrplot(cc, tl.cex = 0.5)

