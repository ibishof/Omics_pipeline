
# Prep data
setwd("C:\\Users\\kqjf682\\Omics_pipelines\\msCleanr\\Erik\\final_data")
data <- read.csv("annotated_MS_peaks-cleaned.csv")
nums <- unlist(lapply(data, is.numeric))
numbers <- data[ , nums]
#numbers <- log(numbers, 2)
clean_data <- cbind(data$Structure, numbers)

# Histogram
names<-names(clean_data)
classes<-sapply(clean_data,class)

par(mfrow=c(2,2))
par(mar=c(5,6,4,2))

for(name in names[classes == 'numeric'])
{
  hist(clean_data[,name],
       main = name,
       breaks=40,
       col="darkmagenta")
}

