

# Your task is to explore this data set and identify any proteins
# that you believe to be changing in abundance through time
# and as a result of the stress induction. Please explain your reasoning.
#
# Please write a brief report which includes code and supporting visualizations
# to address the above questions. Questions do not need to be answered in order
# but clearly indicate what elements of your report address each question.




# Load library
library(dplyr)


############################################################
################### Explore Data ###########################
############################################################


# Bring in data

setwd("C:\\Users\\ibish\\calico")
data <- read.csv("proteomicsSampleData.csv")

# Print unique batch name
unique(data$Batch)
# Four bathes - "HAQ1" "HAQ2" "HAQ3" "HAQ4"

# Print unique time points
unique(data$Time)
# Time points - 20 25 29 33 37 46 50
# Control is missing time point 29

# Print unique Replicates
unique(data$Replicate)
# 1, 2, 3

# Print number of unique peptides
n_distinct(data$Peptide)
# 623 unique peptides

# Group data by ProteinID and summarize the number of unique peptides
protein_single_peptide <- data %>%
  group_by(ProteinID) %>%
  summarise(n_unique_peptide = n_distinct(Peptide)) %>%
  filter(n_unique_peptide == 1)

# Display the ProteinIDs that have only one peptide (One Hit Wonders)
print(protein_peptide_batch$ProteinID)
# sp|P98175-3|RBM10_HUMAN, tr|A4D1B8|A4D1B8_HUMAN 

# Number of unique ProteinIDs
length(unique(data$ProteinID))
# 50 proteins and only 2 are found in a single batch


############################################################
################### Create Reference #######################
############################################################


# Create table that has one row for each peptide
# Use as peptide protein conversion table
# Remove rows with duplicate 'Peptide' values, keeping only the first occurrence
unique_data <- data %>%
  distinct(Peptide, .keep_all = TRUE)

# Create peptide conversion table
peptide_conversion <- select(unique_data, Peptide, ProteinID, Gene)

# Create protein conversion table
protein_conversion <- peptide_conversion  %>% distinct(ProteinID, .keep_all = TRUE)
protein_conversion <- select(protein_conversion, -Peptide)

# Print number of unique Proteins
length(unique(unique_data$Gene))
# 49 unique proteins


############################################################
######## Measure Technical Sources of Variation ############
############################################################


# Calculate variance within each batch for each peptide
batch_variance <- data %>%
  group_by(Peptide, Batch) %>%
  summarise(batch_var = var(Log_Intensity, na.rm = TRUE))

# Calculate variance within each replicate for each peptide
replicate_variance <- data %>%
  group_by(Peptide, Replicate) %>%
  summarise(replicate_var = var(Log_Intensity, na.rm = TRUE))

# Calculate CV for batch variance
mean_batch_var <- mean(batch_variance$batch_var, na.rm = TRUE)
sd_batch_var <- sd(batch_variance$batch_var, na.rm = TRUE)
cv_batch_var <- (sd_batch_var / mean_batch_var)

# Calculate CV for replicate variance
mean_replicate_var <- mean(replicate_variance$replicate_var, na.rm = TRUE)
sd_replicate_var <- sd(replicate_variance$replicate_var, na.rm = TRUE)
cv_replicate_var <- (sd_replicate_var / mean_replicate_var)

# Print out the CVs
print(paste("CV for batch variance: ", cv_batch_var))
print(paste("CV for replicate variance: ", cv_replicate_var))
# CV for batch variance:  2.06159035480846
# CV for replicate variance:  1.59157770830586
# Batch bigger source of variance then replicate


############################################################
############### Hierarchical Models ########################
############################################################


# Hierarchical or mixed-effects models are particularly useful the data is grouped
# observations within the same group are likely correlated.
# Here peptides coming from the same protein can be treated as a group.
# This setup is used to investigate the impact of Time and CellType on Log_Intensity
# while controlling for batch and peptide variability.

# Load libraries
library(lme4)
library(lmerTest)
library(tidyr)
library(stats)

# Initialize empty data frame to store results
results_time_stress <- data.frame(ProteinID = character(), p_value_time_stress = numeric())

# Initialize empty list to store summaries
table_of_summarys <- list()


# Loop through each unique ProteinID
for(protein in unique(data$ProteinID)) {
  
  # Subset data
  data_subset <- subset(data, ProteinID == protein)
  
  # Averaging the replicates for each peptide within each batch, time point, and cell type
  data_subset <- data_subset %>%
    group_by(Peptide, Batch, Time, CellType) %>%
    summarise(Log_Intensity = mean(Log_Intensity, na.rm = TRUE))
  
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


############################################################
####### FDR correct and Make summary table #################
############################################################


# If the Time x Stress interaction is significant less than 0.05 after fdr correction,
# that indicates that the rate of change over time differs between the Stress and Control groups.
# Correct for multiple testing
results_time_stress$p_value_time_stress_adjusted <- p.adjust(results_time_stress$p_value_time_stress, method = "fdr")
results_time_stress$p_value_time_alone_adjusted <- p.adjust(results_time_stress$p_value_time_alone, method = "fdr")

# Rank proteins by adjusted p-value
results_time_stress <- results_time_stress[order(results_time_stress$p_value_time_stress_adjusted), ]

# Add Gene symbol so the Humans can understand
results_time_stress <- merge(protein_conversion, results_time_stress, by = "ProteinID")


