
# Makeing signatures using WGCNA modules

# Create vector with unique colors

library("dplyr")

setwd("~/data/horvath_data2/blood/wgcna")
kMEtable <- read.csv("table_pow8_mergeCutHeight=0.15_kme_bloock800.csv")
colors <- unique(kMEtable$net.colors)
setwd("~/data/horvath_data2")
data <- readRDS("age_balanced_study.rds")
data <- filter(data, age >= 25)

remove(module_summary)
remove(module_avg)
for (i in colors){
  table <- kMEtable[kMEtable$net.colors == i, ]
  module_genes <- table$colnames.proteinGroups.
  module_genes_exp <- dplyr::select(data, module_genes)
  
  # Get rowmean
  module_avg <- rowMeans(module_genes_exp)
  
  if (exists("module_summary")){
    module_summary <- rbind(module_summary, module_avg)
  }
  
  if (isFALSE(exists("module_summary"))){
    module_summary <- module_avg
  }
}

# Add back row and column names
rownames(module_summary) <- colors
colnames(module_summary) <- data$sample

# Transpose and add metadata
module_summary <- as.data.frame(t(module_summary))
module_summary$study <- data$study
module_summary$age <- data$age


# Plot each module and color by study
par(mfrow=c(1,2))
for (i in colors){
  single_color <- dplyr::select(module_summary, age , study, i)
  # single_color <- filter(single_color, study == "GSE42861")
  plot(single_color$age, single_color[,3], main = i, ylab = "Average Methylation", xlab = "age", col = factor(module_summary$study), pch = 20)
}

# Subset interest module for pairs plot
pairs_data <- filter(kMEtable, net.colors == "magenta")
pairs_probes <- pairs_data$colnames.proteinGroups.
pairs_data <- dplyr::select(data, pairs_probes)

pairs(pairs_data, col = factor(module_summary$study), pch = 20)


# Plot 5-10 CpGs per cluster
par(mfrow=c(2,3))
# Number of probes
m <- 6
# Loop though modules and probes
for (i in colors){
  single_module <- dplyr::filter(kMEtable, net.colors ==   i)
  probes <- single_module$colnames.proteinGroups.
  probes <- probes[1:m]
  for (p in probes){
    single_probe <- dplyr::select(data, age, p)
    plot(single_probe$age, single_probe[,2], main = paste(i,p), ylab = "Mean Methylation", xlab = "Age", col = factor(module_summary$study), pch = 20)
  }
  
}

# Calculate PCA for cells color category

# Grab module or modules of interest
setwd("~/data/horvath_data2")
data <- readRDS("age_balanced_study.rds")
data <- filter(data, age >= 25)

setwd("~/data/horvath_data2/blood/wgcna")
kMEtable <- read.csv("table_pow8_mergeCutHeight=0.15_kme_bloock800.csv")
select_colors <- c("magenta")
select_modules <- dplyr::filter(kMEtable, net.colors %in% select_colors)
select_probes <- select_modules$colnames.proteinGroups.
input <- select(data, select_probes)
input[is.na(input)] = 0
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(2,2))
par(mar=c(5,6,4,2))
title <- "PCA labeled by study"



# Make PCA plot with study and print
plot <- autoplot(pca_res,
                 data = data,
                 colour = as.numeric(as.factor(data$study)), # Need '' around column name of legend won't work
                 label = FALSE,
                 label.size = 3,
                 shape = 20,
                 main = title,
)


print(plot)
plot.new()
legend("topleft",   # Coordinates (x also accepts keywords)
       legend = unique(data$study), # Vector with the name of each group
       fill = c(1:length(unique(data$study))),   # Creates boxes in the legend with the specified colors
       # col = par("col"), # Color of lines or symbols
       # border = "black", # Fill box border color
       # lty, lwd,         # Line type and width
       # pch,              # Add pch symbols to legend lines or boxes
       # bty = "o",        # Box type (bty = "n" removes the box)
       # bg = par("bg")    # Background color of the legend
       # box.lwd = par("lwd"), # Legend box line width
       # box.lty = par("lty"), # Legend box line type
       # box.col = par("fg"),  # Legend box line color
       # cex = 1,          # Legend size
       # horiz = FALSE     # Horizontal (TRUE) or vertical (FALSE) legend
       # title = NULL      # Legend title
)


# Calculate PCA for cells color category
input <- age_balanced_top[,-c(1:2)]
input <- as.data.frame(t(input))
input[is.na(input)] = 0
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(2,2))
par(mar=c(5,6,4,2))
title <- "PCA Labeled by Cluster Color"
kMEtable <- kMEtable[order(kMEtable$c.1.ncol.proteinGroups..), , drop = FALSE]

# Make PCA plot with study and print
plot <- autoplot(pca_res,
                 data = input,
                 colour = kMEtable$net.colors, # Need '' around column name of legend won't work
                 label = FALSE,
                 label.size = 10,
                 shape = 19,
                 main = title
)


print(plot)


