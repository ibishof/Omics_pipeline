

time_to_first_reported <- readRDS("~/data/ukbiobank/time_to_first_reported.rds")

# Assuming your dataframe is called time_to_first_reported
# Use colSums to count non-NA values in each column
non_na_counts <- colSums(!is.na(time_to_first_reported))

# Print the result
number_of_diagnosis <- as.data.frame(cbind(colnames(time_to_first_reported), non_na_counts))


library(dplyr)
metadata <- select(ages_diagnosis, eid, age_accel)
time_to_first_reported <- merge(metadata, time_to_first_reported, by = "eid")
age_aceel_per_disease <- list()
for (i in 3:ncol(time_to_first_reported)){
  temp <- time_to_first_reported[,c(2,i)]
  temp <- temp[complete.cases(temp), ]
  age_aceel_per_disease[[i]] <- mean(temp$age_accel)
}

age_aceel_per_disease_table <- as.data.frame(cbind(colnames(time_to_first_reported)[-c(1:2)],
                                                   as.numeric(as.character(unlist(age_aceel_per_disease)))))
age_aceel_per_disease_table$V2 <- as.numeric(as.character(age_aceel_per_disease_table$V2))


age_aceel_per_disease_table <- merge(age_aceel_per_disease_table, number_of_diagnosis, by = "V1")
age_aceel_per_disease_table$non_na_counts <- as.numeric(as.character(age_aceel_per_disease_table$non_na_counts))

age_aceel_per_disease_table <- filter(age_aceel_per_disease_table, non_na_counts > 9 )
