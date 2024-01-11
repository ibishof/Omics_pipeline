library(WGCNA)
library(dplyr)



proteinGroups <- data[,-1]
traitData <- select(data, sample, age )

# Remove rows with too many NaN
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
data[is.nan(data)] <- NA
# View(head(data))
# Remove rows with more than 90% NA
data <- data[which(rowMeans(!is.na(data)) > 0.99), ]

# Remove unwanted metadata
data <- data %>% distinct(sample, .keep_all = TRUE)
sample_id <- data$sample
age <- data$age
data <- data[,-c(1:10)]
probes <- colnames(data)
# View(head(data))

# Clean data
library("DescTools")
dataw = sapply(c(1:ncol(data)), function(x) Winsorize(data[,x], minval = NULL, maxval = NULL, probs = c(0.05, 0.95),  na.rm = TRUE, type = 7) )
data <- as.data.frame(dataw)
colnames(data) <- probes
rownames(data) <- sample_id
hist((sapply(data, var, na.rm = "TRUE")*10000), main = "Var", col = "blue", breaks = 100)
# Remove low variance probes
var0 <- unlist(lapply(data, function(x) 0.001 < var(na.rm = "TRUE", if (is.factor(x)) as.integer(x) else x)))
var0 <- as.logical(as.vector(var0))
data_clean <- data[,var0]

# Fancy transpose
datat <- as.data.frame(t(data_clean))
probes <- colnames(data_clean)
rownames(datat) <- probes
data <- datat


# Choose a set of soft-thresholding powers
# Change proteinGroups to what ever your table is called in this case reduced_compostion
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(proteinGroups, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#you want connecivity below 100 and above 10


###BlockwiseModule function
#netwrok type = signed keeps anti correlate things from being in the same module only time to change is if looking at signaling networks. Things going down and up could be related to a singaling pathway
#Deepslit the higher the number the more modules.
#Pam stageing allows lower corrlated stuff that is still connected on the dentro to be placed in a module (less Gray)
#minModuleSize should be no more than 10% of rows/observations
#mergeCutHeight lower number more modules
allowWGCNAThreads(nThreads = 16)
enableWGCNAThreads(nThreads = 16)
WGCNAnThreads()

net <- blockwiseModules(proteinGroups,power=8,
                        mergeCutHeight=0.15,
                        corType="bicor",
                        #maxPOutliers = 0.05,
                        networkType="signed",
                        pamStage= TRUE,
                        pamRespectsDendro=TRUE, 
                        deepSplit = 2,
                        TOMDenom = "mean",
                        verbose=3,saveTOMs=FALSE,
                        minModuleSize = 3,
                        minKMEtoStay = 0.3,
                        maxBlockSize=800,
                        reassignThreshold = 0.05)


# Plot the dendrogram and the module colors underneath
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = (net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Relating modules to external clinical traits from tutorial
# traitData = read.csv("traits.csv", row.names = 1)
traitData <- traitData %>% distinct(sample, .keep_all = TRUE)
traitData <- filter(traitData, sample %in% sample_id)
traitData <- select(traitData, age)
nGenes = ncol(proteinGroups);
nSamples = nrow(proteinGroups);
MEs = net$MEs
moduleTraitCor = cor(MEs, traitData, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "Age",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#Ordered Modules List of M# and Colors
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))

#Create the meat of the table of the master module kME table:
kMEdat <- signedKME(proteinGroups, net$MEs[,-ncol(net$MEs)], corFnc="bicor")

orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:ncol(proteinGroups)),colnames(proteinGroups),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]

# Add gene symbols
library("FDb.InfiniumMethylation.hg19")

# hm450 <- get450k()
probenames <- kMEtable$`colnames(proteinGroups)`
probes <- hm450[probenames]
tss <- getNearestTSS(probes)
transcript <- getNearestTranscript(probes)
gene_names <- transcript$nearestGeneSymbol
kMEtable$Gene_Name <- gene_names
setwd("~/data/horvath_data2/blood/wgcna")
write.csv(kMEtable,file=paste0("table_pow8_mergeCutHeight=0.15_kme_bloock800.csv"),row.names=FALSE)


############### Boxplots for modules

MEList = moduleEigengenes(t(proteinGroups), colors = net$colors)
MEs = t(MEList$eigengene)


toplot <- as.matrix(MEs)

grouping <- as.data.frame(t(read.csv("groups_ms_notms.csv")))

Group <- as.factor(grouping$V1)


par(mfrow=c(2,2))
par(mar=c(5,6,4,2))

for (i in 1:nrow(toplot)) {
  
  boxplot(toplot[i,]~(Group), names=c("MS","Not_MS"),row=(rownames(MEs))[i],ylab="Eigengene Value",main= (gsub("ME", "", rownames(toplot))[i]),xlab=NULL,las=2, col = (gsub("ME", "", rownames(toplot))[i]) )
  
}
write.csv(toplot, "toplot.csv")






