
# Hot encode drug table

# Split the medications string and create a list of drug names for each participant
drug_list <- strsplit(as.character(drugs$medications), "\\|")

# Get all unique drug names
all_drugs <- unique(unlist(drug_list))
# Remove NA
all_drugs <- na.omit(all_drugs)


# Initialize a matrix to store the one-hot encoded data
hot_encode_matrix <- matrix(0, nrow = nrow(drugs), ncol = length(all_drugs), dimnames = list(NULL, all_drugs))

# Loop through each participant and set the corresponding drug column to 1
for (i in 1:nrow(drugs)) {
  if (!is.na(drugs$medications[i])) {
    hot_encode_matrix[i, drug_list[[i]]] <- 1
  }
}

# Create the hot-encoded dataframe
hot_encode_drug <- as.data.frame(hot_encode_matrix)
hot_encode_drug$eid <- drugs$eid
