## Title: "Proteomics Data Analysis"

## Objective

The task is to explore this data set and identify any proteins that you believe to be changing in abundance through time and as a result of the stress induction. This report includes code and supporting visualizations to address these questions.

# Link with plots and code below:

- [See HTML version with embended plots](https://htmlpreview.github.io/?https://github.com/ibishof/Omics_pipeline/blob/main/Mixed_effects_models/Proteomic_Analysis_Rmarkdown2.html)

```{r setup, include=FALSE}
# Load necessary libraries
library(dplyr)
library(lme4)
library(lmerTest)
library(tidyr)
library(stats)
library(knitr)
library(kableExtra)
```


# Bring in data
```{r}
setwd("C:\\Users\\ibish\\data")
data <- read.csv("proteomicsSampleData.csv")
```


# Explore Data
```{r}
# Print unique batch name
unique(data$Batch)

# Print unique time points
unique(data$Time)

# Print unique Replicates
unique(data$Replicate)

# Print number of unique peptides
n_distinct(data$Peptide)

# Group data by ProteinID and summarize the number of unique peptides
protein_single_peptide <- data %>%
  group_by(ProteinID) %>%
  summarise(n_unique_peptide = n_distinct(Peptide)) %>%
  filter(n_unique_peptide == 1)

# Display the ProteinIDs that have only one peptide (One Hit Wonders)
print(protein_single_peptide$ProteinID)

# Number of unique ProteinIDs
n_distinct(data$ProteinID)

```


# Create Reference
Use as peptide protein conversion table
```{r}
# Remove rows with duplicate 'Peptide' values, keeping only the first occurrence
unique_data <- data %>%
  distinct(Peptide, .keep_all = TRUE)

# Create peptide conversion table
peptide_conversion <- select(unique_data, Peptide, ProteinID, Gene)

# Create protein conversion table
protein_conversion <- peptide_conversion  %>% distinct(ProteinID, .keep_all = TRUE)
protein_conversion <- select(protein_conversion, -Peptide)

# Print number of unique Gene names
n_distinct(unique_data$Gene)
# 49 unique Gene names (Less than ProteinID)
```



# Measure Technical Sources of Variation
CV for batch variance:  0.0604009664033113;
CV for replicate variance:  0.0499613739990236

As expected Batch variance is higher than replicate variance
```{r}
# Calculate SD across batches for each peptide
batch_sd <- data %>%
  group_by(Peptide, Replicate, CellType, Time) %>%
  summarise(batch_sd = sd(Log_Intensity, na.rm = TRUE))

# Calculate variance within each replicate for each peptide
replicate_sd <- data %>%
  group_by(Peptide, Batch, CellType, Time) %>%
  summarise(replicate_sd = sd(Log_Intensity, na.rm = TRUE))

# Calculate means across batches
batch_mean <- data %>%
  group_by(Peptide, Replicate, CellType, Time) %>%
  summarise(batch_mean = mean(Log_Intensity, na.rm = TRUE))

# Calculate means across replicates
replicate_mean <- data %>%
  group_by(Peptide, Batch, CellType, Time) %>%
  summarise(replicate_mean = mean(Log_Intensity, na.rm = TRUE))

cv_batch_var <- mean(batch_sd$batch_sd/batch_mean$batch_mean, na.rm = T)
cv_replicate_var <- mean(replicate_sd$replicate_sd/replicate_mean$replicate_mean, na.rm = T)

# Print out the CVs
print(paste("CV for batch variance: ", cv_batch_var))
print(paste("CV for replicate variance: ", cv_replicate_var))
```

## Running Hierarchical Models on the Data

The following code chunk initializes an empty data frame to store results and loops through each unique ProteinID to fit hierarchical models. Time and treatment (CellType) are treated as fixed effects since I am interested in their relationship to abundance. Batch and Peptide are treated as random effects because I am not interested in them, but I want to account for their possible effects. When Replicates were added as a random effect, the model produced isSingular Fit warnings. These are often due to the random effects structure being too complex. Since Replicates contribute less variance than Batch, I chose to use Batch as the random effect.
```{r}
# Initialize empty data frame to store results
results_time_stress <- data.frame()

# Initialize empty list to store summaries
table_of_summarys <- list()

# Loop through each unique ProteinID for hierarchical models
for(protein in unique(data$ProteinID)) {
   # Subset data
  data_subset <- subset(data, ProteinID == protein)
  
  # Averaging the replicates for each peptide within each batch, time point, and cell type
  data_subset <- data_subset %>%
    group_by(Peptide, Batch, Time, CellType) %>%
    summarise(Log_Intensity = mean(Log_Intensity, na.rm = TRUE), .groups = "keep")
  
  # Check the number of unique batches and peptides for the protein
  # For now ignore proteins with a single peptide or batch
  if(length(unique(data_subset$Batch)) <= 1 || length(unique(data_subset$Peptide)) <= 1) {
    warning(paste("Skipping protein with ID", protein, "because there is only one unique Batch or Peptide."))
    next
  }
  
  # Fit model
  # Batch and Peptide treated as random effect since I am not interested in them but what to account for possible effects
  # When Replicates added as random effect the model did not converge and since they add less variance then batch, I choose batch
    model <- lmer(Log_Intensity ~ Time * CellType + (1|Batch) + (1|Peptide), data = data_subset)

  if (isSingular(model)) {
    warning(paste("Singular fit for protein", protein))
  }
  
  # Extract p-value for interaction term and store summary
  model_summary <- summary(model)
  p_value_time_alone <- model_summary$coefficients[2, "Pr(>|t|)"]
  p_value <- model_summary$coefficients[4, "Pr(>|t|)"]
  
  # Store result and summary
  results_time_stress <- rbind(results_time_stress,
                               data.frame(ProteinID = protein,
                                          p_value_time_alone = p_value_time_alone,
                                          p_value_time_stress = p_value))
  # Save total summaries
  table_of_summarys[[protein]] <- model_summary
}
```


# Correct for multple hypothesis
FDR was used because it provides a balance between controlling for Type I errors (false positives) and not being too conservative so as to miss potentially significant results (Type II errors). Plus commonly used so people are familiar with it.
```{R}
options(scipen=0)
# If the Time x Stress interaction is significant less than 0.05 after fdr correction,
# that indicates that the rate of change over time differs between the Stress and Control groups.
# Correct for multiple testing
results_time_stress$p_value_time_alone_adjusted <- p.adjust(results_time_stress$p_value_time_alone, method = "fdr")
results_time_stress$p_value_time_stress_adjusted <- p.adjust(results_time_stress$p_value_time_stress, method = "fdr")

# Add Gene symbol so the Humans can understand
results_time_stress <- merge(protein_conversion, results_time_stress, by = "ProteinID")

# Rank proteins by adjusted p-value
results_time_stress <- results_time_stress[order(results_time_stress$p_value_time_stress_adjusted), ]

# Make so kable renders as scientific notation. This was super annoying to figure out.
results_time_stress$p_value_time_stress_adjusted <- formatC(results_time_stress$p_value_time_stress_adjusted, format = "e", digits = 2)

```


# Table of p-values
```{R}
kable(results_time_stress, format = "html", table.attr = 'class="table table-striped"') %>%
  column_spec(7, background = "yellow", color = "black")
```


# Correlation Analysis
As an orthogonal approach, I aim to measure the correlation of each peptide with time under both conditions. Subsequently, I will bin the peptides by protein and create plots. Highly significant proteins, such as PDCD11 and SVIL, should exhibit markedly different relationships with time under stressed conditions compared to control conditions.



First, lots of formatting needs to be done to measure correlation. This includes subsetting by condition, summarizing peptide intensities across time points, and handling missing values.


# Subset by Condition
```{R}
# Filter data to subset by condition
control <- filter(data,  CellType == "Control")
stressed <- filter(data, CellType == "Stress")


# Create subsets based on unique time values
for(time in unique(control$Time)) {
  assign(paste0("control_", time), subset(control, Time == time))
}

for(time in unique(stressed$Time)) {
  assign(paste0("stressed_", time), subset(stressed, Time == time))
}
```


# Summarize Peptide Intensities
Ignore Batch and Replicate effects.
```{R}

# Loop through time point and summarize peptide mean intensities
# Group by Peptide and find the mean of each group
# Initialize an empty data frame to store the results
control_means <- data.frame()

# Loop through time points and summarize peptide mean intensities
for (time in unique(control$Time)){
  # Dynamically get the dataframe for this time point
  time_point <- get(paste0("control_", time))
  
  # Group by Peptide and find the mean of each group
  mean_control_time <- time_point %>%
    group_by(Peptide) %>%
    summarize(mean_intensity = mean(Log_Intensity))
  
  # Rename the mean_intensity column to reflect the time point
  colnames(mean_control_time)[which(names(mean_control_time) == "mean_intensity")] <- paste0("mean_control_", time)
  
  # Merge with control_means
  if (nrow(control_means) == 0) {
    control_means <- mean_control_time
  } else {
    control_means <- merge(control_means, mean_control_time, by = "Peptide", all = TRUE)
  }
}

# Now for stressed
# Initialize an empty data frame to store the results
stressed_means <- data.frame()

# Loop through time points and summarize peptide mean intensities
for (time in unique(stressed$Time)){
  # Dynamically get the dataframe for this time point
  time_point <- get(paste0("stressed_", time))
  
  # Calculate the mean intensities for this time point
  mean_stressed_time <- time_point %>%
    group_by(Peptide) %>%
    summarize(mean_intensity = mean(Log_Intensity))
  
  # Rename the mean_intensity column to reflect the time point
  colnames(mean_stressed_time)[which(names(mean_stressed_time) == "mean_intensity")] <- paste0("mean_stressed_", time)
  
  # Merge with stressed_means
  if (nrow(stressed_means) == 0) {
    stressed_means <- mean_stressed_time
  } else {
    stressed_means <- merge(stressed_means, mean_stressed_time, by = "Peptide", all = TRUE)
  }
}


```


# Explore Peptide Completeness
Control has more missing data then stress
```{R}

# Are there missing values?
hist(rowSums(is.na(control_means[,-1])), main = "Number of NA's per row control", xlab = "Number of NA", col = "blue")
# Number of NA per column control
colSums(is.na(control_means[,-1]))

# Plot with histogram
hist(rowSums(is.na(stressed_means[,-1])), main = "Number of NA's per row stressed", xlab = "Number of NA", col = "blue")
# Number of NA per column stress
colSums(is.na(stressed_means[,-1]))
```


# Remove peptides with half or more of the values being NA
```{R}
# Count the number of NA values for each row, excluding the first column
na_counts <- apply(control_means[ , -1], 1, function(x) sum(is.na(x)))

# Identify rows where half or more of the values are NA (ignoring the ID column)
rows_to_remove <- na_counts >= (ncol(control_means) - 1) / 2

# Remove those rows
clean_control_means <- control_means[!rows_to_remove, ]

# Now for stress condition
na_counts <- apply(stressed_means[ , -1], 1, function(x) sum(is.na(x)))
rows_to_remove <- na_counts > (ncol(stressed_means) - 1) / 2
clean_stressed_means <- stressed_means[!rows_to_remove, ]

```


# Find correlation between all features and age
```{R}

# Find correlation between all features and age
time_points_control <- unique(control$Time)

# Calculate correlation between control peptides and time
c_control = sapply(c(1:nrow(clean_control_means)), function(x) cor(time_points_control, unlist(clean_control_means[x,-1]),
                                            method = "pearson", use = "complete.obs") )

# Get p-values for correlation, this account for number of data points and strength of correlation in one value
p_control = sapply(c(1:nrow(clean_control_means)), function(x) cor.test(time_points_control, unlist(clean_control_means[x,-1]), method = "pearson", use = "complete.obs") )

# Transform the list to a data frame for easier manipulation and readability.
p_control <- as.data.frame(t(as.data.frame(p_control)))
p_control <- unlist(p_control$p.value)

# FDR correct
# Adjust for multiple test. Since a large number of peptides compared to time points. To avoid Type II errors FDR was choosen.
library(stats)
p_control_fdr <- p.adjust(as.numeric(p_control), method = "fdr")


# Calculate correlation between stressed peptides and time
time_points_stress <- unique(stressed$Time)
c_stressed = sapply(c(1:nrow(clean_stressed_means)), function(x) cor(time_points_stress, unlist(clean_stressed_means[x,-1]),
                                                               method = "pearson", use = "complete.obs") )

# Get p-values for correlation, this account for number of data points and strength of correlation in one value
p_stressed = sapply(c(1:nrow(clean_stressed_means)), function(x) cor.test(time_points_stress, unlist(clean_stressed_means[x,-1]), method = "pearson", use = "complete.obs") )

# Transform the list to a data frame for easier manipulation and readability.
p_stressed <- as.data.frame(t(as.data.frame(p_stressed)))
p_stressed <- unlist(p_stressed$p.value)

# FDR correct
p_stressed_fdr <- p.adjust(as.numeric(p_stressed), method = "fdr")

```


# Summarize correlation and p-values into nice readable tables
```{R}
# Initialize a data frame to store correlation metrics for peptides.
# Each metric (e.g., correlation coefficient, p-value, FDR-adjusted p-value, etc.) is placed in a separate column.
peptide_control_corr_table <- data.frame(
  Peptide          = clean_control_means$Peptide,
  Correlation      = as.numeric(c_control),
  PValue           = p_control,
  FDR_Adjusted_P   = p_control_fdr,
  NegLog10_FDR_P   = -log10(unlist(p_control_fdr)),
  Num_NAs          = apply(clean_control_means[, -1], 1, function(x) sum(is.na(x)))
)

peptide_control_corr_table <- merge(peptide_conversion, peptide_control_corr_table, by = "Peptide")

# For stressed cells
peptide_stressed_corr_table <- data.frame(
  Peptide          = clean_stressed_means$Peptide,
  Correlation      = as.numeric(c_stressed),
  PValue           = p_stressed,
  FDR_Adjusted_P   = p_stressed_fdr,
  NegLog10_FDR_P   = -log10(unlist(p_stressed_fdr)),
  Num_NAs          = apply(clean_stressed_means[, -1], 1, function(x) sum(is.na(x)))
)


peptide_stressed_corr_table <- merge(peptide_conversion, peptide_stressed_corr_table, by = "Peptide")

```


# Hand Check highly correlates peptides
```{R}

# Filter the row and drop the Peptide column
single_peptide <- clean_stressed_means %>% filter(Peptide == "K.IYEDGDDDMKR.T") %>% select(-Peptide)

# Convert the row to a numeric vector
single_peptide <- as.numeric(unlist(single_peptide))

plot(time_points_stress, single_peptide, pch = 19, col = "blue", main = "K.IYEDGDDDMKR.T")

```


# Boxplots of peptide correlations
Each peptide is plotted as a single data point. The "spread" in correlation for each protein becomes easy to visualized.
```{R}
# Sort peptide_stressed_corr_table by Gene A-Z
peptide_stressed_corr_table <- peptide_stressed_corr_table[order(peptide_stressed_corr_table$Gene), ]

# Sort peptide_control_corr_table by Gene A-Z
peptide_control_corr_table <- peptide_control_corr_table[order(peptide_control_corr_table$Gene), ]


par(mfrow = c(2, 4))
for ( protein in unique(peptide_stressed_corr_table$ProteinID)){
  
  # Subset data
  data_subset_stressed <- subset(peptide_stressed_corr_table, ProteinID == protein)
  data_subset_control <- subset(peptide_control_corr_table, ProteinID == protein)
  
  # Create boxplot
  boxplot(data_subset_stressed$Correlation, ylim = c(-1,1), main = paste0(unique(data_subset_stressed$Gene), " stressed"), col = "red", ylab = "Corr")
  boxplot(data_subset_control$Correlation, ylim = c(-1,1), main = paste0(unique(data_subset_control$Gene), " control"), col = "blue", ylab = "Corr")
}
```


# Top 2 Most Significantly Changed and Least Significantly Changed Proteins Due to Time and Stress
Here, one can clearly see that the correlation patterns for peptides belonging to proteins PDCD11 and SVIL are very different in stressed (red) and control (blue) conditions. In contrast, PSMD5 and RB1CC1 look very similar under both control and stressed conditions. The correlations for their peptides are also highly variable, moving both up and down, suggesting that these proteins do not exhibit a consistent linear relationship to time in this experiment.

```{R}

gene_names <- c("PDCD11", "SVIL", "PSMD5", "RB1CC1")

# Filter peptide_stressed_corr_table by gene_name
peptide_stressed_corr_table <- filter(peptide_stressed_corr_table, Gene %in% gene_names)

# Filter peptide_control_corr_table by gene_name
peptide_control_corr_table <- filter(peptide_control_corr_table, Gene %in% gene_names)

par(mfrow = c(2, 4))
for (protein in gene_names){
  
  # Subset data
  data_subset_stressed <- subset(peptide_stressed_corr_table, Gene == protein)
  data_subset_control <- subset(peptide_control_corr_table, Gene == protein)
  
  # Create boxplot
  boxplot(data_subset_stressed$Correlation, ylim = c(-1,1), main = paste0(unique(data_subset_stressed$Gene), " stressed"), col = "red", ylab = "Corr")
  boxplot(data_subset_control$Correlation, ylim = c(-1,1), main = paste0(unique(data_subset_control$Gene), " control"), col = "blue", ylab = "Corr")
}

```
