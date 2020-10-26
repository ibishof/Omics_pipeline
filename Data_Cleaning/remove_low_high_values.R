
setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\NIH\\casey")
proteinGroups<-read.csv("substract.csv", row.names = 1)
mat <- as.matrix(proteinGroups)


mat[mat > 0.0] <- 0


write.csv(mat, "negative.csv")

write.csv(mat, "postive.csv")
