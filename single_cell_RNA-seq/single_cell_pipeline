
# Seurat/SingleR single cell pipeline

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(sctransform)

# For annotation
library(SingleR)
library(celldex)





  # For local files:
  setwd("C:\\Users\\ibish\\Data\\rna_seq\\single_cell\\GSE249669")
  expression_matrix <- ReadMtx(
    mtx = "GSE249669_matrix.mtx", 
    features = "GSE249669_genes.tsv", # Often called genes
    cells = "GSE249669_barcodes.tsv"
  )
 
  # For remote files:
  expression_matrix <- ReadMtx(mtx = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE249669&format=file&file=GSE249669%5Fmatrix%2Emtx%2Egz",
                               cells = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE249669&format=file&file=GSE249669%5Fbarcodes%2Etsv%2Egz",
                               features = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE249669&format=file&file=GSE249669%5Fgenes%2Etsv%2Egz")
  seurat_object <- CreateSeuratObject(counts = data)



# Initialize the Seurat object with the raw (non-normalized data).
# Use the count matrix to create a Seurat object. 
# The object serves as a container that contains both data (like the count matrix)
# and analysis (like PCA, or clustering results) for a single-cell dataset.
# For more information, check out our [Seurat object interaction vignette], or our GitHub Wiki.
# Gene must be found in at least 3 cells
# Cell must have 200 features or more to be included
pbmc <- CreateSeuratObject(counts = expression_matrix,
                           project = "pbmc3k",
                           min.cells = 3,
                           min.features = 200)
pbmc
# 24485 Gens across 109927 Cells

# Examine the count matrix
counts <- GetAssayData(pbmc, layer = 'counts')
counts[c("SNRNP70", "TCL1A", "MS4A1", "CD19"), 1:50]


# Time for QC #
# The number of unique genes detected in each cell.
# Low-quality cells or empty droplets will often have very few genes
# Cell doublets or multiplets may exhibit an aberrantly high gene count
# Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
# The percentage of reads that map to the mitochondrial genome
# Low-quality / dying cells often exhibit extensive mitochondrial contamination
# We calculate mitochondrial QC metrics with the PercentageFeatureSet() function, which calculates the percentage of counts originating from a set of features
# We use the set of all genes starting with MT- as a set of mitochondrial genes

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
MT_pecentage <- as.data.frame(pbmc[["percent.mt"]])


# Visualize QC metrics as a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships plot(x,y), but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells that have unique feature counts over 2,500 or less than 200
# A very low count of unique features might indicate a cell with poor quality RNA
# Cells with an unusually high number of unique gene counts might indicate doublets or multiplets. 
# These are instances where two or more cells have been captured as a single "cell"

# Filter cells that have >5% mitochondrial counts
# Cells with a high percentage of mitochondrial gene expression are often indicative of stressed or dying cells.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 5)        
pbmc
# 24485 features across 2610 cells


# Normalizing the data #

# First step is to normalize the gene expression levels in each cell by the total expression in that cell, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result.
# This step transforms the data into something akin to counts per ten thousand (CPTT)
# Gene expression is  right-skewed, a few genes are very highly expressed while most genes have low expression. 
# The log transformation reduces the skewness, making the data more normally distributed. 
# This normalization makes downstream analysis techniques, like clustering or principal component analysis, 
# more effective and meaningful.

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection #

# Identify genes that show significant variation in expression levels across different cells.
# Since this is count data mean expression corresponds to variance
# Identify gene with higher variance then thier mean expression would correspond too
# These genes likly represent genes that define cell types, states, or responses to environmental stimuli
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data #
# Get z-score
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# Remove unwanted sources of variation
pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt")

# Perform linear dimensional reduction #
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Major contributors to PCs
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# PCA plot
DimPlot(pbmc, reduction = "pca", cols = "blue") + NoLegend()

# Setting cells to a number plots the ‘extreme’ cells on both ends of the spectrum,
# which dramatically speeds plotting for large datasets.
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# Plot top 15 PCs
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# What percent of Variance in captured in each PC?
# Repeat downstream analyses with a different number of PCs (10, 15, or even 50!). Check if results differ dramatically
ElbowPlot(pbmc)

# Cluster the cells #
# K-nearest neighbor (KNN) graph, with edges drawn between cells with similar feature expression patterns,
# and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.(Louvain method)
# Calculate Jaccard similarity using FindNeighbors() function

pbmc <- FindNeighbors(pbmc, dims = 1:15)

# Between 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells.
# Optimal resolution often increases for larger datasets.
pbmc <- FindClusters(pbmc, 
                     resolution = 0.5) # Increased values leading to a greater number of clusters.

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# UMAP/tSNE
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pbmc <- RunUMAP(pbmc, dims = 1:15) # dims = number of PCs used
DimPlot(pbmc, reduction = "umap")

# Save analysis
setwd("C:\\Users\\ibish\\Data\\rna_seq\\single_cell\\GSE249669")
saveRDS(pbmc, file = "pbmc.rds")
pbmc <- readRDS("pbmc.rds")

# Find Cluster Biomarkers #

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Top biomarkers for each cluster
top_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  arrange(p_val) %>%
  slice_head(n = 5)

# Get AUC values for different genes
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Visualize cell type biomarkers
biomarkers <- top_markers$gene
VlnPlot(pbmc, features = biomarkers)


# Plot raw counts
VlnPlot(pbmc, features = biomarkers, slot = "counts", log = TRUE)

# Examine cell clusters via known biomarkers
# Native t-cell VlnPlpt
VlnPlot(pbmc, c("IL7R", "CCR7", "SELL"))

# TH2 t-cells
VlnPlot(pbmc, c("GATA3", "IL4", "IL2","IL13"))

# TH1 t-cells
VlnPlot(pbmc, c("IFNG", "TBX21", "IL5","IL12RB1"))

# Interferon-related immune active t-cells
VlnPlot(pbmc, c("OAS1", "OAS2", "OAS3", "EIF2AK2"))

# Memory CD4+
VlnPlot(pbmc, c("IL7R", "S100A4"))

# CD14+ Mono VlnPlpt
VlnPlot(pbmc, c("CD14", "LYZ"))

# FCGR3A+ Mono
VlnPlot(pbmc, c("FCGR3A", "MS4A7"))

# CD8+T
VlnPlot(pbmc, "CD8A")

# B-cells VlnPlpt
VlnPlot(pbmc, c("MS4A1"))

# NK cells
VlnPlot(pbmc, c("GNLY", "NKG7"))

# DC cells
VlnPlot(pbmc, c("FCER1A", "CST3"))

# Platelet VlnPlpt
VlnPlot(pbmc, c("PPBP"))


# DoHeatmap() generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Annotation
# Tutorial here: https://www.youtube.com/watch?v=7RuPGaWcY0Y

ref <- celldex::HumanPrimaryCellAtlasData()
# View labels
View(as.data.frame(colData(ref)))
counts <- GetAssayData(pbmc, layer = 'counts')

# Perform annotation
singleR_results <- SingleR(test = counts, ref = ref, labels = ref$label.main)
singleR_results

# Add annotation to Seurat object
pbmc$singleR.labels <- singleR_results$labels[match(rownames(pbmc@meta.data), rownames(singleR_results))]
DimPlot(pbmc, reduction = 'umap', group.by = 'singleR.labels')

# Add SingleR annotations to Seurat object metadata
pbmc <- AddMetaData(pbmc, metadata = singleR_results$labels, col.name = "singleR_pred")

# Diagnostics
# Columns are cells and rows are labels
plotScoreHeatmap(singleR_results)

# See how many cell in each category was pruned
plotDeltaDistribution(singleR_results)
summary(is.na(singleR_results$pruned.labels))

# Examine the expression of the marker genes for each label in the dataset
all.markers <- metadata(singleR_results)$de.genes
count$labels <- singleR_results$labels

# Compare results to unsupervised clustering
tab <- table(Assigned = singleR_results$labels, Clusters = pbmc$seurat_clusters)
tab
# make heat map, add 10 to aviod NAs
#
pheatmap::pheatmap(log10(tab+10), color = colorRampPalette(c('white', 'blue'))(10))





# Assigning cell type identity to clusters #
new.cluster.ids <- c("Activated CD8+T", # 0
                     "Naive T-cell", # 1
                     "CD14+ Mono", # 2
                     "NK cells", # 3
                     "T-cells TH1", # 4
                     "T-cells active", # 5
                     "CD14+ Mono", # 6
                     "Interferon T-cells", # 7
                     "T-cells TH2", # 8
                     "Trained Monocytes", # 9
                     "B cell", # 10
                     "DC", # 11 
                     "Platelet") # 12

levels(pbmc) <- as.character(c(1:length(levels(pbmc))-1))
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
ggsave(filename = "../umap.jpg", height = 7, width = 12, plot = plot, quality = 50)


# find all markers distinguishing Native t-cell from clusters 0 and 3
cluster1.markers <- FindMarkers(pbmc, ident.1 = "t-cells 2", ident.2 = c("t-cells TH1", "t-cell TH2", "Native t-cell"))
head(cluster1.markers, n = 10)


cluster0.markers <- FindMarkers(pbmc, ident.1 = "CD8+T pure", ident.2 = c("t-cells TH1", "CD8 T"))
head(cluster0.markers, n = 10)
VlnPlot(pbmc, "NME1")
