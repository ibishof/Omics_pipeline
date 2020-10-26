


install.packages("dplyr")
install.packages("d3heatmap")
install.packages("gplots")
install.packages("ggrepel")
install.packages("devtools")
install.packages("ggbiplot")



library(dplyr)
library(RColorBrewer)
library(d3heatmap)
library(ggplot2)
library(gplots)
library(ggrepel)
library(devtools)
library(ggbiplot)

setwd("C:\\Users\\bishofij\\Desktop")

exp1 = read.csv("proteomicSeminar.csv")
#Checkout data set
head(exp1)
dim(exp1)

newdata <- exp1[c(4,20:22)]
head (newdata)
str(newdata)
# get column names
colnames(newdata)

# Rename column by name
names(newdata)[names(newdata) == "Abundance.Ratio...F1..Medium.....F1..Light."] <- "R1"
names(newdata)[names(newdata) == "Abundance.Ratio...F2..Medium.....F2..Light."] <- "R2"
names(newdata)[names(newdata) == "Abundance.Ratio...F3..Medium.....F3..Light."] <- "R3"
head(newdata)

#Do I have NA or any other string factors in my date frame?
table(is.na(newdata))

sum(is.na(newdata$R1))
sum(is.na(newdata$R2))
sum(is.na(newdata$R3))

#Where my NA are?
missingvalues = newdata [which(is.na(newdata$R2)),]
missingvalues

#Dropping variables 
dropping_newdata <- subset(newdata, R1>= 0.02 & R2>= 0.02 & R3>= 0.02 & R1<= 90 & R2<= 90 & R3<= 90 ) 
dim(dropping_newdata)
#Do I have NA or any other string factors in my date frame?
table(is.na(dropping_newdata))

sum(is.na(dropping_newdata$R1))
sum(is.na(dropping_newdata$R2))
sum(is.na(dropping_newdata$R3))
dropping_newdata

write.csv(dropping_newdata, "dropping_newdata0.02and99.csv")
#how many observation have been dropped?
dim (newdata)
#How many where lost
dim(dropping_newdata)

#heatmap 
datamatrix=(as.matrix(dropping_newdata[, -1]))
heatmap.2(datamatrix,col=brewer.pal(11, "RdYlBu"),scale="row", trace="none")

#Can not due PCA need ggbiplot
si.log2.pca <- prcomp(dropping_newdata[,c(2:4)], center = TRUE,scale. = TRUE)
summary(si.log2.pca)
ggbiplot(si.log2.pca) +geom_point(colour = "deepskyblue4", size = 0.05) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                panel.background = element_blank(), axis.line = element_line(colour = "black"))

#ttest for volcano plot
data_for_analysis<-dropping_newdata[,c(1,2:4)]
t_test<-function(x)p_val=t.test(x,mu=1,na.rm=T)$p.value
columns<-data_for_analysis[,2:4]
ttest<-as.vector(apply(columns,1,t_test))
print(ttest)

Finalttest = cbind(dropping_newdata, ttest)
colnames(Finalttest)
head(Finalttest)
write.csv(ttest, "Finalttest.csv")
p.adjust = p.adjust(Finalttest$ttest, method = "fdr")

finalFDR = cbind(Finalttest, p.adjust)
head(finalFDR)

names(finalFDR )[names(finalFDR ) == "ttest"] <- "Raw.pvalue"
names(finalFDR )[names(finalFDR ) == "p.adjust"] <- "FDR"
head(finalFDR )

log2newdata = log2(dropping_newdata[2:4])
head(log2newdata)

Average.R = apply(log2newdata[, 1:3], 1, mean)
head(Average.R)

log2final = cbind(finalFDR , Average.R)
colnames(log2final)

volcanoplot = log2final[c(1,7,6)]

head (volcanoplot)
write.csv(volcanoplot, "volcanoplot.csv")

#Graph volcano plot
volcano = ggplot(data = volcanoplot, aes(x = Average.R, y = -log10(FDR)))

volcano + geom_point(alpha=0.3, position = position_jitter(), size = 0.8, colour= "grey") + coord_cartesian(ylim = c(0, 2.5)) +geom_hline(yintercept = 1.3, color = "hotpink", size =0.2) + geom_vline(xintercept = -2.5, color = "hotpink", size =0.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) + geom_text_repel(size=2.5, data=subset(volcanoplot, Average.R < -3.5 & FDR < 0.05),segment.size  = 0.1,
                                                                                                          segment.color = "hotpink", aes(Average.R,-log10(FDR), label = Accession), check_overlap = TRUE) +  xlab("Log2FoldChange") + ylab("-Log10(FDR)")
