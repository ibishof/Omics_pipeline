##############################################################################
## One-Step GO-ELITE WITH USER PARAMETERS - by Eric Dammer, Divya Nandakumar
##  - performs GO-Elite v1.2.5 with Fisher Exact Test for enrichment p<0.05
##    and 5 minimum genes per ontology
##
## Nicholas Seyfried Lab Bioinformatics - for the lab - 02/11/2019 update
##  - Change parameters between the lines below
##############################################################################
### GO-Elite is a python package to perform ontology analysis on gene sets. 
### The python script for GO-Elite can be downloaded from http://www.genmapp.org/go_elite/
### Alternatively, there is also a GUI which can be downlaoded from the same website.
###
### Custom databases can be downloaded from http://software.broadinstitute.org/gsea/msigdb/
### including TRANSFAC and curated pathways collected from various sources, e.g. the C2 DB.
### GO-Elite requires python 2.7 be installed, and FET with command-line requires
### that a bugfix be applied; copy GO_Elite.py included to GO-Elite program subfolder:
### GOeliteFolder/GO-Elite_v.1.2.5-Py/
##############################################################################

rm(list = ls())
dev.off()
options(stringsAsFactors=FALSE)

######################## EDIT THESE VARIABLES (USER PARAMETERS) #########################################################
fileName <- "master_module_kME8.csv"    #Input list, often comes from WGCNA pipeline
            #INPUT is a CSV FILE - in the filePath folder.(The variable below)
            #Can be formatted as Kme table from WGCNA pipeline, or
            #can be a CSV of columns, one symbol or UniqueID (Symbol|...) list per column, with the LIST NAMEs in row 1
            #in this case, the longest list is used as background for GO-Elite.
            #  For simple columnwise list input, DON'T FORGET TO PUT THE APPROPRIATE BACKGROUND LIST IN, OR RESULTS WILL BE UNRELIABLE.
filePath <- "C:/Users/bishofij/Proteomics_Pipeline/go-elite_r_projects/"
            #Folder that will contains the input file specified above, and where your output file will be created. 
            #for system call from R to python to work, avoid folders with spaces; 
            #but in Windows, spaces are handled by the script, making them passable.

outFilename <- "test6"
            #Where you want the output to go. SUBFOLDER WITH THIS NAME WILL BE CREATED AUTOMATICALLY
GOeliteFolder <- "C:/Users/bishofij/Proteomics_Pipeline/"
            #has subfolder and python script for GO-Elite v1.2.5 (script should be edited per authors' instructions to Divya Nandakumar to correct bug for Fisher Exact Test
            #Remove from end of file path GO-Elite_v.1.2.5-Win64/GO-Elite_v.1.2.5-Py
            #for system call from R to python to work, avoid folders with spaces; 
            #but in Windows, spaces are handled by the script, making them passable.

maxBarsPerOntology=5
            #Ontologies per ontology type, used for generating the PDF report; does not limit GO-Elite output
speciesCode="Hs"
            #Hs for homo sapiens, Dm for fly; Mm, mouse; Rn, rat... (must have database downloaded via command line)
            #if you use the GUI for GO-Elite, create a copy of the folders for command line, and delete databases downloaded via GUI.
downloadDB=FALSE
            #If TRUE, the database files for speciesCode species will be downloaded "(be patient)" from Ensembl (v62 preferred)
pythonPath <- "C:/Python27/"
            #python.exe for python v2.7 is here.

panelDimensions=c(3,2)
            #rows, columns for each page of PDF output
color=c("darkseagreen3","lightsteelblue1","lightpink4")
            #colors respectively for ontologyTypes:
            #"Biological Process","Molecular Function","Cellular Component"
            #must be valid R colors
modulesInMemory=FALSE
            #uses cleanDat, net, and kMEdat from pipeline already in memory
ANOVAgroups=FALSE
            #if true, modulesInMemory ignored. Volcano pipeline code should already have been run!
############ ADVANCED OPTION ####################################################################################
customDBcmd=paste0("--customSet ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/Databases/EnsMart62Plus/C2/ --dataToAnalyze all ")
            #set to "" if you have no custom database.
customPanelDimensions=c(2,2)
customReportMaxBars=20
######################## REST OF THE CODE IS AUTOMATIC ##########################################################


## Clean out spaces and escaped backslashes from folder paths (folder names with spaces should not be used on non-windows systems with this script)
filePath=paste0(paste( sapply(do.call(c,strsplit(filePath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
pythonPath=paste0(paste( sapply(do.call(c,strsplit(pythonPath,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")
GOeliteFolder=paste0(paste( sapply(do.call(c,strsplit(GOeliteFolder,"[/\\]")),function(x) { if (grepl(" ",x)) { gsub(x,paste0(substr(gsub(" ","",x),1,6),"~1"),x) } else { x } } ),collapse="/"),"/")


## The input files for GO-Elite are text files with the gene list as the 1st column, a symbol identified (gene symbol, uniprot etc) as the 2nd column
## Different accepted inputs are given in the tutorial
## Commonly used symbols - Gene Symbol - Sy (example of input file below)
### GeneSymbol		SystemCode (Symbol format)
###	  GFAP		Sy
###	  APOE		Sy
## All input files are placed in one folder

## The background file is prepared similarly and is placed in a separate folder
## The initial part of the code prepares files for GO-Elite. This can be skipped if the files are being made manually as described above.
## The second part of the code runs GO-ELite either from R (using the system command) or can be run using the terminal (in mac)
## The second part requires GO-Elite to be installed and path to the GO-Elite installation site indicated following python
## The 3rd part of the code plots the results from the GO-Elite results folder. When using the GUI the 1st 2 parts can be skipped and only the 3rd part can be used for plotting

##-------------------------------##
## Preparing files for GO-Elite ##
## Takes in the module assignment file as input with 1st column having gene names, 2nd column having color assignments followed by kME values


dir.create(file.path(filePath, outFilename))

##1a. GO Elite of the significant up and down (p<0.05) proteins in the current cleanDat
### GO Elite analysis for ANOVA-defined categories ###
if (ANOVAgroups) {
  sigThresh=0.05;

  #colnames(ANOVAout)
  #********************
  ##This code relies on a pre-existing data frame ANOVAout already in the environment, processed through volcano output for selected pairwise comparisons of interest.
  #numComp=6 #number of pairwise comparisons for ANOVA+Tukey p value columns, which are followed by the same-order log2(mean difference) columns
  #********************

  ANOVAout$Symbol <- do.call("rbind", strsplit(as.character(rownames(ANOVAout)), "[|]"))[, 1]

  DEXlistsForGO<-list()
  iter=0
  for (i in testIndexMasterList) {
    iter=iter+1;
    j=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".down")
    k=paste0(gsub(" ",".",comparisonIDs$Comparison[iter]),".up")
    if (length(intersect(i,flip))==1) {
      #flipped sign (all diffs >0 for down)
      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]>0)]
      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]<0)]
    } else {
      #do not flip sign (all diffs <0 for down)
      DEXlistsForGO[[j]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]<0)]
      DEXlistsForGO[[k]]<-ANOVAout$Symbol[which(ANOVAout[,i]<sigThresh & ANOVAout[,i+numComp]>0)]
    }
  }

  #write lists to GOElite input files, and also the background file
  for (i in names(DEXlistsForGO)) { 
    dfGO<-data.frame(GeneSymbol=DEXlistsForGO[[i]],SystemCode=rep("Sy",length(DEXlistsForGO[[i]])))
    write.table(unique(dfGO),file=paste(filePath,outFilename,"/",i,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
  }
  #write background
  background <- unique(ANOVAout$Symbol)
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  dir.create(file.path(paste0(filePath,outFilename),"background"))
  write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
  nModules=length(names(DEXlistsForGO))
  WGCNAinput=FALSE

} else { #NOT creating lists from ANOVA/volcano up & down groups


##1b. GO Elite of the WGCNA modules from specified input file or in the current cleanDat, net, and kME table
### GO Elite analysis for WGCNA-defined modules ###

##use data structures in memory if modulesInMemory=TRUE; otherwise read csv following the template which is written in an earlier R session and edited in excel to produce module membership table saved as .csv:
if (modulesInMemory) {
  modulesData <- cbind(rownames(cleanDat),net$colors,kMEdat)
  WGCNAinput=TRUE
} else {
  modulesData <- read.csv(paste(filePath,fileName,sep=""),header=TRUE, sep=",");
  # check if this is a WGCNA modules/kME table or simple list input
  if(length(na.omit(match("net.colors",colnames(modulesData))))>0) {
    WGCNAinput=TRUE
    # Remove first and last column if it contains sort information ("Original Order" and "color|kMEin")
    modulesData <- modulesData[,-c(1,ncol(modulesData))];
  } else {
    WGCNAinput=FALSE
  }
}

if (WGCNAinput) {
  library(WGCNA) #for labels2colors

  # Include column with Symbol (if it is gene symbol, if not use appropriate code as given in GO-Elite manual)
  modulesData$SystemCode <- rep("Sy",nrow(modulesData)) 

  # Assign Names of First columns, in case they are non standard
  colnames(modulesData)[1]<-"Unique.ID" #This should have Symbol|UniprotID
  colnames(modulesData)[2]<-"net.colors" #This should have colors

  #Split out symbols from UniprotIDs, keep symbols in column 1
  rownames(modulesData)<-modulesData$Unique.ID
  modulesData$Unique.ID<-do.call("rbind",strsplit(as.character(modulesData$Unique.ID), "[|]"))[,1]

  ## Creating background file for GO Elite analysis
  background <- unique(modulesData[,"Unique.ID"])
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  dir.create(file.path(paste0(filePath,outFilename),"background"))
  write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

  # Separate into independent module txt files for analysis by GO-Elite (CREATE INPUT FILES)
  greySubtractor=if(length(which(modulesData$net.colors=="grey"))>0) { 1 } else { 0 } #remove grey from count of modules
  nModules <- length(unique(modulesData$net.colors))-greySubtractor
  moduleColors <- uniquemodcolors <- labels2colors(c(1:nModules)) 
  for (i in 1:length(moduleColors)) {
    moduleName <- moduleColors[i]
    ind <- which(colnames(modulesData) == gsub("kMEME","kME",paste("kME",moduleName,sep="")))
    moduleInfo <- modulesData[modulesData$net.colors == gsub("ME","",moduleName), c(1,ncol(modulesData),ind)]
    colnames(moduleInfo) <- c("GeneSymbol","SystemCode","kME")
    if (moduleName == "blue" | moduleName == "brown" | moduleName == "green" | moduleName == "cyan") { write.table(moduleInfo,file=paste(filePath,outFilename,"/",moduleName,"_2_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    } else {
      write.table(unique(moduleInfo),file=paste(filePath,outFilename,"/",moduleName,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
    }
  }
} else { #input is not WGCNA kME table format

##1c. GO Elite of the WGCNA modules from specified input file, which must be in a column-wise list format, and including longest such list as background.
  # We process the input file as simple lists by column in the CSV (largest list used as background)

  #reread the file to a list of gene symbol (or UniqueID) lists
  modulesData <- as.list(read.csv(paste(filePath,fileName, sep=""),sep=",", stringsAsFactors=FALSE,header=T)) 

  nModules <- length(names(modulesData))
  for (a in 1:nModules) {
    modulesData[[a]] <- unique(modulesData[[a]][modulesData[[a]] != ""])
    modulesData[[a]] <- modulesData[[a]][!is.na(modulesData[[a]])]
    modulesData[[a]] <- do.call("rbind",strsplit(as.character(modulesData[[a]]), "[|]"))[,1]
  }
  ## Creating background file for GO Elite analysis
  background <- modulesData[order(sapply(modulesData,length),decreasing=TRUE)][[1]]
  background <- unique(background)
  background <- cbind(background,rep("Sy",length=length(background)))
  colnames(background) <- c("GeneSymbol","SystemCode")
  dir.create(file.path(paste0(filePath,outFilename),"background"))
  write.table(background,paste0(filePath,outFilename,"/background/background.txt"),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

  # Separate Symbol Lists into independent module txt files for analysis by GO-Elite (CREATE INPUT FILES)
  modulesData[[ names(modulesData[order(sapply(modulesData,length),decreasing=TRUE)])[1] ]] <- NULL
  nModules = nModules -1 #no background
  listNames <- uniquemodcolors <- names(modulesData)
  for (i in listNames) {
    listName <- i
    listInfo <- cbind(modulesData[[listName]],rep("Sy",length=length(modulesData[[listName]])))
    colnames(listInfo) <- c("GeneSymbol","SystemCode")
    write.table(unique(listInfo),file=paste(filePath,outFilename,"/",listName,".txt",sep=""),row.names=FALSE,col.names=TRUE,sep="\t", quote=FALSE)
  }
} #end if (WGCNAinput)
} #end else for if (ANOVAgroups)


##2. GO Elite Python Call
####----------------------- GO-ELite Analysis in Command Prompt/ Terminal or from R ------------------------------------#####

# Input and Denominator Files were prepared in the specified format (3 columns - genelist, format (Sy- for gene symbol) and kME)
# Each module is a separate txt file. Denominator file contains background list whcih is all identified proteins in this case
# All input files are placed in input or geneInfo folder and background file in denominator folder

if (downloadDB) {
  cat(paste0("Downloading Ensembl v62 database for species code ",speciesCode,". (Be patient.)\n\n"))
  system( gsub("'",'\"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --update Official --species ",speciesCode," --mod Ensembl --version EnsMart62Plus")) )
}

commandLine=gsub("'",'"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --species ",speciesCode," --mod Ensembl --permutations 'FisherExactTest' --method 'z-score' --zscore 1.96 --pval 0.05 --num 5 --input ",filePath,outFilename,"/ --denom ",filePath,outFilename,"/background/ ",customDBcmd,"--output ",filePath,outFilename,"/"))
cat(paste0("NOW RUNNING THE FOLLOWING COMMAND:\n\n", commandLine,"\n\n(Estimated time for ", nModules, " lists to complete: ",round((30*nModules)/60,1)," minutes)\n Start time: ",Sys.time(),"\n"))
system( gsub("'",'\"',paste0(pythonPath,"python ",GOeliteFolder,"GO-Elite_v.1.2.5-Py/GO_Elite.py --species ",speciesCode," --mod Ensembl --permutations 'FisherExactTest' --method 'z-score' --zscore 1.96 --pval 0.05 --num 5 --input ",filePath,outFilename,"/ --denom ",filePath,outFilename,"/background/ ",customDBcmd,"--output ",filePath,outFilename,"/")) )


##3. Output Report of Z-Score Barplots, processing all GO-Elite output files
############################# ----------------------Plotting for modules ------------------------#######################
######## this script plots the top 3 ontologies for biological process, mol function and cell component for each module


##color scheme for ontology type key/legend (can be changed in user parameters, editing the "color" vector)
ontologyTypes=c("Biological Process","Molecular Function","Cellular Component")

if(ANOVAgroups) {
  xlabels <- names(DEXlistsForGO)
  xlabels.frame <- data.frame(Colors=rep(NA,length(xlabels)),Labels=xlabels)
  uniquemodcolors <- names(DEXlistsForGO) #not set above
} else {
  xlabels <- uniquemodcolors #labels2colors(c(1:nModules))
  xlabels1 <- paste("M",seq(1:nModules),sep="")
  xlabels.frame <- as.data.frame(data.frame(Colors=xlabels,Labels=paste0(xlabels1," ",xlabels)))
}

setwd(paste0(filePath,outFilename,"/"))
pdf(paste0("GO-Elite_",outFilename,".pdf"),height=10,width=8)
op <- par(mfrow=panelDimensions,oma=c(0,0,3,0))
frame()
legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)
#frame()
GOEliteOUTfileTrailer<-if(ANOVAgroups | !WGCNAinput) { c("-GO_z-score_elite.txt"); } else { c("_Module-GO_z-score_elite.txt"); }
summary <- list()
for(i in c(1:(length(uniquemodcolors)))){
	thismod=uniquemodcolors[i]	
	if (file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep="")) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
		if(file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep="")) == T) {
			tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep=""),sep="\t")
		} else {
		next
		} 
	} else {
	tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep=""),sep="\t") #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
	}
	if (length(tmp[,2]) == 0) next
	tmp = tmp[,c(2,3,9,10,11)] ## Select GO-terms,GO-Type, Z-score,pValues and gene Lists
	tmp1 = tmp[order(tmp$Z.Score,decreasing=T),]
	tmp2 = tmp1[order(tmp1$Ontology.Type,decreasing=T),] #was tmp2
	tmp3 = tmp2[tmp2$Ontology.Type == "biological_process",][c(1:maxBarsPerOntology),]
	tmp3 = rbind(tmp3,tmp2[tmp2$Ontology.Type == "molecular_function",][c(1:maxBarsPerOntology),] )
	tmp3 = rbind(tmp3,tmp2[tmp2$Ontology.Type == "cellular_component",][c(1:maxBarsPerOntology),] )
	tmp3 <- na.omit(tmp3)
#	tmp3 <- tmp3[order(tmp3$Z.Score,decreasing=T),] #added this row, if you want to mix ontology types and sort by Z.Score only
	tmp3 <- tmp3[rev(rownames(tmp3)),]

	summary[[i]] <- tmp3
	
	### To color bars by mol function, cell component or biological process
	for (j in 1:nrow(tmp3)){
		if (tmp3$Ontology.Type[j] == "molecular_function"){
			tmp3$color[j] <- color[2]
		} else if (tmp3$Ontology.Type[j] == "cellular_component"){
			tmp3$color[j] <- color[3]
		} else if (tmp3$Ontology.Type[j] == "biological_process"){
			tmp3$color[j] <- color[1]
		}
	# tmp3$color[j] <- uniquemodcolors[i] #module color for all bars, instead of different colors by ontology type
	} 

	if (tmp3$Z.Score == F) next
	par(mar=c(4,15,4,3))
	xlim <- c(0,1.1*max(tmp3$Z.Score))	
	moduleTitle <- xlabels.frame[i,"Labels"]
	xh <- barplot(tmp3$Z.Score,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=0.9,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
	abline(v=1.96,col="red", cex.axis = 0.5)
	axis(2, at=xh, labels = tmp3$Ontology.Name, tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
}

par(op) # Leaves the last plot
dev.off()


if(!customDBcmd=="") {
library(stringr)
pdf(paste0("GO-Elite_",outFilename,"-CUSTOM_db.pdf"),height=10,width=8)
op <- par(mfrow=customPanelDimensions,oma=c(0,0,3,0))
frame()
legend(x="topleft",legend = ontologyTypes, fill=color, title=" ",cex=2,horiz=F,xpd=T)
legend(x="topleft",legend = c(" "," "," "), title="Ontology Types",cex=2.5,horiz=F,xpd=T, bty='n', title.adj=1.4)

GOEliteOUTfileTrailer<-if(ANOVAgroups | !WGCNAinput) { c("-UserSuppliedAssociations_z-score_elite.txt"); } else { c("_Module-UserSuppliedAssociations_z-score_elite.txt"); }
summary <- list()
for(i in c(1:(length(uniquemodcolors)))){
	thismod=uniquemodcolors[i]	
	if (file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep="")) == F) { #**note "_Module" needs to be removed from fileTrailer if setting uniquemodcolors not by module color
		if(file.exists(paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep="")) == T) {
			tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,"_2",GOEliteOUTfileTrailer,sep=""),sep="\t")
		} else {
		next
		} 
	} else {
	tmp=read.csv(file=paste(filePath,outFilename,"/GO-Elite_results/CompleteResults/ORA_pruned/",thismod,GOEliteOUTfileTrailer,sep=""),sep="\t")
	}
	if (length(tmp[,2]) == 0) next
	tmp = tmp[,c(1,1,7,8,12)] ## Select GO-terms,GO-Type, Z-score,pValues and gene Lists
	tmp1 = tmp[order(tmp$Z.Score,decreasing=T),]
	tmp2 = tmp1
	tmp3 = tmp2
	tmp3 <- na.omit(tmp3)
	tmp3 <- tmp3[order(tmp3$Z.Score,decreasing=T),][c(1:customReportMaxBars),] #added this row, if you want to mix ontology types and sort by Z.Score only
	tmp3 <- tmp3[rev(rownames(tmp3)),]

	summary[[i]] <- tmp3
	
	### To color bars by mol function, cell component or biological process
	for (j in 1:nrow(tmp3)){
		tmp3$color[j] <- color[1]
		if(WGCNAinput) { tmp3$color[j] <- uniquemodcolors[i] } #module color for all bars, instead of different colors by ontology type
	} 

#	if (tmp3$Z.Score == F) next
	if (is.na(max(tmp3$Z.Score))) tmp3<-na.omit(tmp3)

	par(mar=c(4,15,4,3))
	xlim <- c(0,1.1*max(tmp3$Z.Score))	
	moduleTitle <- xlabels.frame[i,"Labels"]
	xh <- barplot(tmp3$Z.Score,horiz = TRUE,width =0.85,las=1,main=moduleTitle, xlim=xlim,col=tmp3$color,cex.axis=0.7,xlab="Z Score",cex.lab=0.9,cex.main=0.95,ylim=c(0,nrow(tmp3)+0.8))
	abline(v=1.96,col="red", cex.axis = 0.5)
	axis(2, at=xh, labels = str_to_title(gsub("_"," ",tmp3$Gene.Set.Name)), tick=FALSE, las =2, line =-0.5, cex.axis = 0.7) 
}

par(op) # Leaves the last plot
dev.off()
} #end if(!customDBcmd=="")

