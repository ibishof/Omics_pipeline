
# Prep data
setwd("C:\\Users\\kqjf682\\Omics_pipelines\\CKD\\A003_1000")
data <- read.csv("Area__noQC_blank_no_ignore.csv")
nums <- unlist(lapply(data, is.numeric))
numbers <- data[ , nums]

numbers <- log(numbers, 2)
is.na(numbers) <- do.call(cbind,lapply(numbers, is.infinite))
numbers[is.na(numbers)] <- 0

clean_data <- cbind(data$Metabolite.name, numbers)
clean_numbers <- numbers
numbers <- t(numbers)


library(matrixStats)
library(dplyr)
# Calculate the mean and median of each sample
no_ID <- numbers
data_mean <- rowMeans(no_ID)
data_median <- rowMedians(as.matrix(no_ID, na.rm = TRUE))
data_zero <- apply(no_ID == 0, 1, sum) 

# Create bar graph
par(mfrow=c(1,3))
par(mar=c(5,6,4,2))
barplot(data_zero, main="Zeros",
        xlab="Samples")
barplot(data_mean, main="Mean",
        xlab="Samples")
barplot(data_median, main="Median",
        xlab="Samples")



# Histogram
names<-names(clean_numbers)
classes<-sapply(clean_numbers,class)

par(mfrow=c(3,3))
par(mar=c(5,6,4,2))

for(name in names[classes == 'numeric'])
{
  hist(t(dv),
       main = "Correlation Histogram",
       breaks=100,
       col="darkmagenta")
}

par(mfrow=c(4,3))
par(mar=c(5,6,4,2))
for (col in 1:ncol(clean_numbers)) {
  hist(clean_numbers[,col],
       breaks = 40,
       col = "darkmagenta")
}


library(ggplot2)
library(ggfortify)
# Remove zero variance columns
input<-numbers
input <- input[ , which(apply(input, 2, var) != 0)]

# Calculate PCA
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
#title <- paste0("Test ",as.numeric(ncol(training)))

# Make PCA plot and print
plot <- autoplot(pca_res,
                 data = input,
                 #colour = 'Diagnosis',
                 #label= input$V1,
                 shape = FALSE,
                 #frame = TRUE
                 main = "PCA")
print(plot)



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

