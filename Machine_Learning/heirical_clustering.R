

clusdata <- read.csv("combat.csv", row.names = 1)

# Makes simple dentrogram
clusdata <- t(clusdata)
d <- dist(clusdata, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
plot(hc1, cex = 0.6, hang = -1, las= .025, cex = .1, edge.width = 0.1 )


# Fancy looking dentrogram
library(dendextend)
library(colorspace)

dend <-  as.dendrogram(hc1)

# Bring in file containing col names and groups
car_type <- read.csv("levels.csv")
car_type <- car_type$Cohort
car_type <- factor(car_type)
n_car_types <- length(unique(car_type))
# Turns factors into R colors
cols_2 <- colorspace::rainbow_hcl(n_car_types, c = 70, l  = 50)
# Matches colors to car type
col_car_type <- cols_2[car_type]

# color labels by car company:
labels_colors(dend) <- col_car_type[order.dendrogram(dend)]
# color branches based on cutting the tree into 4 clusters:
#dend <- color_branches(dend, k = 2)

### plots
par(mar = c(12,4,1,1))
# Set label size
dend <- set(dend, "labels_cex", 0.35)
plot(dend, horiz = TRUE)
#colored_bars(cbind(k234[,3:1], col_car_type), dend1, rowLabels = c(paste0("k = ", 4:2), "Car Type"))

