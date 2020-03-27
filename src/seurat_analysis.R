library(dplyr)
library(Seurat)
library(Matrix)


# Load the  dataset
data <- Read10X(data.dir = "./data/filtered_gene_bc_matrices/mm10/")
# Initialize the Seurat object with the raw (non-normalized data).
ob <- CreateSeuratObject(counts = data, project = "DO", min.cells = 3, min.features = 200)


#percent MT
ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^mt-")

# Visualize QC metrics as a violin plot
VlnPlot(ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(ob, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

ob <- subset(ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5800 & percent.mt < 10)


ob <- NormalizeData(ob, normalization.method = "LogNormalize", scale.factor = 10000)

#identify highly variable features (high cell-to-cell variation)
ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(ob), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(ob)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

#Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
all.genes <- rownames(ob)
ob <- ScaleData(ob, features = all.genes)

ob <- RunPCA(ob,features = VariableFeatures(object = ob))

# Examine and visualize PCA results a few different ways
print(ob[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(ob, dims = 1:2, reduction = "pca")

DimPlot(ob, reduction = "pca")

#Determine the ‘dimensionality’ of the dataset

#jackstraw: randomly permute 1% of data, and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. 
#We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.

ob <- JackStraw(ob, num.replicate = 100)
ob <- ScoreJackStraw(ob, dims = 1:20)

JackStrawPlot(ob, dims = 1:20)

ElbowPlot(ob)


##
ob <- FindNeighbors(ob, dims = 1:12)
ob <- FindClusters(ob, resolution = 0.5)
##

ob <- RunUMAP(ob, dims = 1:12)

DimPlot(ob, reduction = "umap", label=T)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

FeaturePlot(ob, features = c("Rasd1"))

#save(ob, file = "./results/Rdata/seurat_ob.Rdata")

DimPlot(ob, reduction = "umap", label=TRUE,split.by = "seurat_clusters")

cluster5.markers <- FindMarkers(ob, ident.1 = 0, min.pct = 0.25)

ob.markers <- FindAllMarkers(ob, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ob.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

#markers distinguishing cluster
cluster5.markers <- FindMarkers(ob, ident.1 = c(5), ident.2 = c(0,1,2,3,4,6,7,8,9,10), min.pct = 0.25)
head(cluster5.markers, n = 30)

# FeaturePlot(ob, features = c("Qsox1"))
# 
# 
# FeaturePlot(ob, features = c("Tnfsf11", "Tnfrsf11a", "Tnfrsf11b"))
# 
# 
# FeaturePlot(ob, features = c("Sp7", "Ibsp", "Bglap", "Col1a1"))
# 
# FeaturePlot(ob, features = c("Sost", "Dmp1", "Phex", "Mepe"))
# 
# FeaturePlot(ob, features = c("Ctsk", "Ocstamp", "Dcstamp"))
# FeaturePlot(ob, features = c("Adcy9"))
# FeaturePlot(ob, features = c("Rasd1"))
