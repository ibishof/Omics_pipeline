# Omic pipeline
R scripts for proteomic and metabolic data analysis

**Machine Learning:**
A series of machine learning algorithms used for regression and classification. 
- **New Tree Based Pipelines:**
These pipelines use either Random forest, Extra Treees, or XGboost to both select features and build models. The pipelines build a model (Random forest, Extra Trees, XGboost), calculates accuratcy, feature importance, plots seperation via PCA, and recursively selects features. 

**MaxQuant Bash Scripts:**
Scripts for running MaxQuant on a HPC in a Linux enviroment.

**QC:**
This contains a QC pipeline that can be used to idenify outliers and bad runs/samples in a dataset. Calculates and graphs number of zeros, median, and mean instenisty across samples. Histrograms of each feature, PCA plot and correlation matrix. 

**WGCNA:**
Weighted correlation network analysis (WGCNA) is used for finding clusters (modules) of highly correlated genes/protiens. The central idea is that proteins that are correlated have some type of biological relatedness. Part of the pipeline, line 112, produces a table that is then used for GO-elite.

**GO-elite:**
The objective of GO-elite is to identify a  set of biological Ontology terms or pathways to describe a particular set of genes/proteins. It answeres the basic question who are these proteins and what do they do. The code here is an R based wrapper for GO-elite that is used to visulize the results.

**Random Forest:**
The objective of this pipeline is to build a model that can predict diagnosis based of relative protein abundance. The pipeline also calculates feature importance. This information can be used to find biomarker candidates.

**New Tree Based Pipelines:**
These pipelines use either Random forest, Extra Treees, or XGboost to both select features and build models. PCA is then used to visualize the seperation of classes. For regression analysis the true vs predicted values are plotted.

**R_scripts:**
This contains short scripts rangeing from data cleaning, QC, model building, and visualization of results.

Together this pipelines can inform the researcher what group of proteins are related to disease. What pathways these proteins are in and what proteins are best used as biomarkers.
