


# Calculate PCA of Features
library(ggfortify)
input <- small_data[,-c(1:3)]
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
#title <- paste0("Test ",as.numeric(ncol(training)))



# Make PCA plot and print
library(ggplot2)
library(ggfortify)
plot <- autoplot(pca_res,
                 data = small_data,
                 colour = 'brief_icd10',
                 label = FALSE,
                 label.size = 4,
                 shape = 19,
                 #frame = TRUE,
                 main = "PCA")
print(plot)

