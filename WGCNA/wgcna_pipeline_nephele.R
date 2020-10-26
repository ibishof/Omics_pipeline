#Installing WGCNA package not necissary if previosuly installed

source("https://bioconductor.org/biocLite.R")
biocLite("BiocUpgrade")
source("https://bioconductor.org/biocLite.R")
biocLite()
install.packages("BiocManager") 
BiocManager::install("WGCNA") 

#Load WGCNA

library(WGCNA)

#Load in file
options(stringsAsFactors = FALSE);
setwd("C:\\Users\\bishofij\\Proteomics_Pipeline\\WGCNA\\data_sets")
proteinGroups <-read.csv("proteinGroups_filtered.csv", row.names = 1)

# The point of this block of code is to choose a set of soft-thresholding powers
# Change proteinGroups to what ever your table is called
powers = c(c(1:10), seq(from = 12, to=20, by=2))
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
# This line corresponds to using an R^2 cut-off of h
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
net <- blockwiseModules(t(proteinGroups),power=8,
                        mergeCutHeight=0.2, #detectCutHeight=0.9999, 
                        corType="bicor",networkType="signed",pamStage= TRUE,
                        pamRespectsDendro=TRUE, 
                        deepSplit = 2,
                        TOMDenom = "mean",
                        verbose=3,saveTOMs=FALSE,
                        minModuleSize = 20,
                        minKMEtoStay = 0.3,
                        maxBlockSize=10000,reassignThreshold = 0.05)


# Plot the dendrogram and the module colors underneath
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = (net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Relating modules to external clinical traits from tutorial
#Trait data has to be provided
traitData = read.csv("traits_IP.csv", row.names = 1)
nGenes = ncol(proteinGroups);
nSamples = nrow(proteinGroups);
MEs = net$MEs
moduleTraitCor = cor(MEs, t(traitData), use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = rownames(traitData),
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
kMEdat <- signedKME(t(proteinGroups), net$MEs[,-ncol(net$MEs)], corFnc="bicor")

#Table produced here is used for the Go Elite pipeline
orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(proteinGroups)),rownames(proteinGroups),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
write.table(kMEtable,file=paste0("master_module_kME8.txt"),sep="\t",row.names=FALSE)


############### Boxplots for modules

MEList = moduleEigengenes(t(proteinGroups), colors = net$colors)
MEs = t(MEList$eigengene)


toplot <- as.matrix(MEs)

#The grouping is the experminet groups such as Control and Disease
grouping <- as.data.frame(t(read.csv("group_IP.csv")))

Group <- as.factor(grouping$V1)


par(mfrow=c(2,2))
par(mar=c(5,6,4,2))
#names are hard coded in currently
for (i in 1:nrow(toplot)) {
  
  boxplot(toplot[i,]~(Group), names=c("CTL","WT", "dLC1", "dLC2", "dLC1/2"),row=(rownames(MEs))[i],ylab="Eigengene Value",main= (gsub("ME", "", rownames(toplot))[i]),xlab=NULL,las=2, col = (gsub("ME", "", rownames(toplot))[i]) )
  
}







