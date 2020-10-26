
# Crate vector with unique colors
colors <- unique(kMEtable$`net$colors`)

# Bring in table to be subseted
data <- read.csv("znorm.csv")

# Make a table for each module
# Make a list for features in each module
for (i in colors){
  table <- kMEtable[kMEtable$`net$colors` == i, ]
  list <- table$`rownames(proteinGroups)`
  
  
  assign(paste0(i, "_names"), list)
  assign(i, table)
  
  random_forest_input <-data %>% 
    select(Diagnosis,list)
  assign(paste0(i, "_random_forest_input"), random_forest_input)
}


