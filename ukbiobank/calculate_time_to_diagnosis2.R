
library(dplyr)

diagnosis_dates2 <- read.csv("~/data/ukbiobank/diagnosis_dates2.csv")
dates <- read.csv("~/data/ukbiobank/dates.csv")

dates <- dates %>% arrange(eid)
diagnosis_dates2 <- diagnosis_dates2 %>% arrange(eid)

# Save column name for later
names <- colnames(diagnosis_dates2)


# Convert date strings to Date objects
dates$baseline_date <- as.Date(dates$Date.of.attending.assessment.centre...Instance.0, format = "%m/%d/%Y")
diagnosis_dates2 <- lapply(diagnosis_dates2[,-1], function(x) as.Date(x, format = "%m/%d/%Y", origin = "1970-01-01"))

# Calculate time difference in years for each column
time_to_diagnosis <- mapply(function(diagnosis_date) {
  as.numeric(difftime(diagnosis_date, dates$baseline_date, units = "days")) / 365.25
}, diagnosis_dates2)

# Combine results into a table
time_to_diagnosis <- as.data.frame(time_to_diagnosis)

# Add eid column to time_to_diagnosis
time_to_diagnosis$eid <- dates$eid

# Reorder columns with eid as the first column
time_to_diagnosis <- time_to_diagnosis[, c("eid", names(diagnosis_dates2))]

# Rename columns
colnames(time_to_diagnosis) <- names

# View the resulting table
setwd("~/data/ukbiobank")
write.csv(time_to_diagnosis, "time_to_diagnosis.csv")

