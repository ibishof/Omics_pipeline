

# Calculate age when blood was taken
setwd('~/data/ukbiobank')
dates <- read.csv("dates.csv")

# Convert month names to numeric values
dates$month_num <- match(tolower(dates$Month.of.birth), tolower(month.name))

# Convert birth and assay dates to Date objects
dates$birth_date <- as.Date(paste(dates$Month.of.birth, "15", dates$Year.of.birth, sep = " "), format = "%B %d %Y")
dates$Date.of.attending.assessment.centre...Instance.0 <- as.Date(dates$Date.of.attending.assessment.centre...Instance.0, format = "%m/%d/%Y")

# Calculate time difference in years
dates$age_at_blood_draw <- as.numeric(dates$Date.of.attending.assessment.centre...Instance.0 - dates$birth_date) / 365.25


combo <- merge(dates, age, by = "eid")
write.csv(combo, "age_at_blood_drawn.csv")