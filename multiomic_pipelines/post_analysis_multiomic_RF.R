# Post analysis of the features selected by the multiomic RF pipeine
# Find correlation between features and AUC, cluster by indication, and perform PCA
# Find the correlation between features and response to drug

# First select features, many option are provided
# Option One
# Hand Selected related features
features_list <- c("gene1",
                   "gene2",
                   "gene3",
                   "gene4"
)

# Option two use top features from a model iteration
# Select top features
setwd("/home/ibishof/analysis/")
test_cor <- read.csv("results_summarised.csv")
features <- read.csv("58_variable_importance_calc.csv")
features_list <- features$Features

# Option three use top feature from "small data" set
# The "small data" contains the top features selected across multiple models
setwd("/home/ibishof/analysis/PARG/eurofin/ic50_log10")
small_data <- read.csv("small_data.csv")
features_list <- colnames(small_data[,-c(1:2)])
features_list <- gsub("\\-", "\\.", features_list)
features_list <- gsub("\\|", "\\.", features_list)

# Option four use top feature directly from small data set where NA's are kept
setwd("/home/ibishof/analysis/PARG/eurofin/ic50_log10")
smaller_data <- read.csv("small_data_all_true.csv")
smaller_data <- smaller_data[!duplicated(smaller_data[,'DepMap_ID']),]


# Next Read in data
# Option one
# Read in all data except transcripts
setwd("/home/ibishof/ccle/ready")
presmall_data <- readRDS("all_data_minus_transcripts.rds")

# Prep data Remove repeat rows
n_occur <- data.frame(table(presmall_data$DepMap_ID))
single_occurance <- n_occur[n_occur$Freq == 1,]
single_occurance <- as.vector(single_occurance$Var1)
presmall_data <- filter(presmall_data, DepMap_ID %in% single_occurance)
names <- as.vector(colnames(presmall_data))
names <- gsub("gd_gd_", "gd_", names)
colnames(presmall_data) <- names

iRNA <- readRDS("D2_combined_iRNA.rds")
iRNA$CCLE_Name <- rownames(iRNA)

setwd("/home/ibishof/ccle")
master_table <- read.csv("master_name_table2.csv")
iRNA <- merge(master_table, iRNA, by = "CCLE_Name")
col_names <- colnames(iRNA)
colnames(iRNA) <- paste("iRNA", col_names, sep = "_")
names(iRNA)[3] <- "DepMap_ID"
presmall_data <- merge(presmall_data, iRNA, by = "DepMap_ID", all = TRUE)

# Option two
# Read in specific data set
setwd("/home/ibishof/ccle/ready")
smaller_data <- readRDS("expdata_ready.rds")

# Read in transcripts
setwd("/home/ibishof/ccle/ready")
transcripts <- readRDS("transcripts_ready.rds")
# Select features that start with transcripts
top_transcripts <- features_list[(grep("transcripts_", features_list))]
col_names <- colnames(transcripts)
colnames(transcripts) <- paste("transcripts", col_names, sep = "_")
names(transcripts)[1] <- "DepMap_ID"
# Select transcripts needed for analysis
transcripts_small <- select(transcripts, "DepMap_ID", top_transcripts)
# Merge it all together
presmall_data <- merge(presmall_data, transcripts_small, by = "DepMap_ID", all = TRUE)

# Select top features
smaller_data <- presmall_data %>%
  select(DepMap_ID, features_list)


# Read in metadata many option available
# Option one
setwd("/home/ibishof/ccle")
metadata <- read.csv("metadata2_PARG.csv")

# Option two
# If using Niraparib
metadata_niraparib <-niraparib.metdata(metadata_niraparib)
# calculate difference between PARP and PARG
PARP_minus_PARG <- metadata_niraparib$AUC - metadata_niraparib$AUC_IDB001516
metadata_niraparib$PARP_minus_PARG <- PARP_minus_PARG
metadata_niraparib_big_diff <- filter(metadata_niraparib, PARP_minus_PARG > .1 | PARP_minus_PARG < -.1 )
metadata <- metadata_niraparib_big_diff
metadata <- select(metadata, DepMap_ID, PARP_minus_PARG)

# Option three  add eurofin metadata
# Pull in metadata
setwd("/home/ibishof/ccle")
metadata <- read.csv("Eurofins_US034-0013363_Project_27Oct2021_preped.csv")
metadata <- metadata %>% 
  rename(
    CCLE_Name = CCLE_name
  )
metadata <- filter(metadata, compoundId == "IDB001516" )
metadata$ic50_log10 <- log10(metadata$IC50)

# Once metadata is read in
################ Variable ################
Metadata_variables <- c("DepMap_ID", "ic50_log10")
metadata <- select(metadata, Metadata_variables)

# Add metadata
smaller_data_meta <- merge(metadata, smaller_data, by = "DepMap_ID")
vars_remove <- c("DepMap_ID")
smaller_data_meta <- select(smaller_data_meta, -vars_remove)
smaller_data_meta <- filter(smaller_data_meta, !is.na(ic50_log10))
data <- smaller_data_meta

# Remove NA
data[is.na(data)] = 0
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data[is.nan(data)] <- 0


# Begining of analysis
# Find correlation between all features and drug
drug <- data$ic50_log10

c = sapply(c(2:ncol(data)), function(x) cor(drug, data[,x],
                                            method = "pearson", use = "complete.obs") )
# Merge correlation values vector with Gene names
drug_corr_table <- as.data.frame(cbind(colnames(data[,-1]),c))
drug_corr_table$c <- as.numeric(as.character(drug_corr_table$c))
num <- as.numeric(as.character(drug_corr_table$c))
setwd("/home/ibishof/analysis/PARG/eurofin/ic50_log10")
write.csv(drug_corr_table, "expression_feature_corr.csv")

# Create histogram
# dev.off()
hist(num, col = "blue", breaks = 30)


# Calculate PCA of Features
library(ggfortify)
library(caret)
input <- data[,-c(1:2)]
input[is.na(input)] = 0
input <- as.data.frame(t(input))
input <- input[,-nearZeroVar(input)]
input[ , which(apply(input, 2, var) != 0)]
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
#title <- paste0("Test ",as.numeric(ncol(training)))

# Make PCA plot and print
plot <- autoplot(pca_res,
                 data = input,
                 # colour = 'AUC_IDB001516',
                 label = FALSE,
                 label.size = 2.5,
                 shape = TRUE,
                 #frame = TRUE
                 main = "PCA")
print(plot)



# Calculate PCA for cells color by senstivity
input <- data[,-c(1:2)]
input[is.na(input)] = 0
input <- input[,nearZeroVar(input)]
input <- input[ , which(apply(input, 2, var) != 0)]
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
#title <- paste0("Test ",as.numeric(ncol(training)))

# Make PCA plot and print
plot <- autoplot(pca_res,
                 data = data,
                 colour = 'ic50_log10',
                 label = FALSE,
                 label.size = 3,
                 shape = 19,
                 #frame = TRUE,
                 main = "PCA")
print(plot)


# Calculate PCA for cells and color by indication
# Add metadata

vars_retain <- c("DepMap_ID", "primary_disease", "PARP_minus_PARG")
smaller_data_meta <- select(metadata_niraparib_big_diff, vars_retain)
data <- merge(smaller_data_meta, smaller_data, by = "DepMap_ID")


# Remove NA
data[is.na(data)] = 0
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data[is.nan(data)] <- 0

# Remove low frquency indications
freq <- data %>%
  group_by(primary_disease) %>%
  summarize(length(primary_disease))
freq <- filter(freq, `length(primary_disease)` > 12)
freq_indications <- freq$primary_disease
data <- filter(data, primary_disease %in% freq_indications)

# PCA of cell lines and colored by indications
input <- data[,-c(1:2)]
input[is.na(input)] = 0
pca_res <- prcomp(input, scale. = TRUE)
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
#title <- paste0("Test ",as.numeric(ncol(training)))

# Make PCA plot and print
plot <- autoplot(pca_res,
                 data = data,
                 colour = 'primary_disease',
                 label = FALSE,
                 label.size = 4,
                 shape = 19,
                 #frame = TRUE,
                 main = "PCA")
print(plot)


# Build table with correlation of each feature to the AUC per indication
vars <- colnames(data[,-c(1:2)])
remove(one_gene)
remove(ddcr)
for (i in vars){
  one_gene <- data %>%
    group_by(primary_disease) %>%
    summarize(cor(PARP_minus_PARG, get(i)))
  
  if (exists("ddcr")){
    ddcr <- merge(ddcr, one_gene, by = "primary_disease")
  }
  if (isFALSE(exists("ddcr"))){
    ddcr <- one_gene
    
  }
}
colnames(ddcr) <- c("primary_disease", vars)

# Create heat map of features
# Using DDCR file
library(corrplot)
library(ggplot2)
com <-  select(ddcr, -primary_disease)
com[is.na(com)] = 0
com <- com[ , which(apply(com, 2, var) != 0)]


# Now create a correlation matrix with your community composition data using the command 'cor':
cc = cor(com, method = "spearman")
###if you have missing values in your data add the 'use' parameter. Check '?cor' for the cor command help page
# I usually use Spearman correlation because I'm not overly concerned that my relationships fit a linear model,
#and Spearman captures all types of positive or negative relationships (i.e. exponential, logarithmic)


# Plot correlation
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
corrplot(cc, tl.cex = 0.2)


col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = cc, col = col, symm = TRUE, cexRow = .7, cexCol = .2)

# CLustered features in order
heat_object <- heatmap(x = cc, col = col, symm = TRUE, cexRow = 1, cexCol = .5)
old_order <- colnames(com)
new_order <- old_order[heat_object$colInd]
new_order
new_order <- as.data.frame(new_order)



# Heatmap For indications
# Now create a correlation matrix with your community composition data using the command 'cor':
com <- t(com)
colnames(com) <- ddcr$primary_disease
cc = cor(com, method = "spearman")
###if you have missing values in your data add the 'use' parameter. Check '?cor' for the cor command help page
# I usually use Spearman correlation because I'm not overly concerned that my relationships fit a linear model,
#and Spearman captures all types of positive or negative relationships (i.e. exponential, logarithmic)


# Plot correlation
par(mfrow=c(1,1))
par(mar=c(5,6,4,2))
corrplot(cc, tl.cex = 0.2)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = cc, col = col, symm = TRUE, cexRow = 1, cexCol = .5)
heat_object <- heatmap(x = cc, col = col, symm = TRUE, cexRow = 1, cexCol = .5)
# CLustered features in order
old_order <- colnames(com)
new_order <- old_order[heat_object$colInd]
new_order
# Write and save everything
setwd("/home/ibishof/analysis/PARP/PARP_minus_PARG_big_diff/small_data")
write.csv(drug_corr_table, "top_expression_feature_corr_to_auc.csv")
write.csv(ddcr, "corr_expression_by_indication.csv")
