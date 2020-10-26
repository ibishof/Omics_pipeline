
# PCA
library(ggfortify)
training <- read.csv("test.csv", check.names=FALSE)
input<- training[,-1]

# Remove zero variance columns
input$Diagnosis = as.numeric(as.factor(input$Diagnosis))
input <- input[ , which(apply(input, 2, var) != 0)]
pca_res <- prcomp(input, scale. = TRUE)

# Plot
autoplot(pca_res, data = training, colour = 'Diagnosis')



