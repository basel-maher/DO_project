library(dplyr)
library(Seurat)
library(Matrix)
library(data.table)

set.seed(8675309)




s.genes <- cc.genes$s.genes

g2m.genes <- cc.genes$g2m.genes

gene_and_ensemb <- fread("./data/filtered_gene_bc_matrices/mm10//genes.tsv", header=FALSE)
s.genes_ensembl <- gene_and_ensemb[toupper(gene_and_ensemb$V2) %in% s.genes,]
s.genes_ensembl <- s.genes_ensembl$V2
g2m.genes_ensembl <- gene_and_ensemb[toupper(gene_and_ensemb$V2) %in% g2m.genes,]
g2m.genes_ensembl <- g2m.genes_ensembl$V2
#Check cell cycle effects:
############################################################################# Initialize the Seurat object with the raw (non-normalized data)
# Load the  dataset
data <- Read10X(data.dir = "./data/filtered_gene_bc_matrices/mm10/")
# Initialize the Seurat object with the raw (non-normalized data).
ob <- CreateSeuratObject(counts = data, project = "DO", min.cells = 3, min.features = 200)

ob[["percent.mt"]] = PercentageFeatureSet(object = ob, pattern = "^mt-")
ob <- subset(x = ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5800 & percent.mt < 10)

#Scaling the data:
ob <- SCTransform(ob)


ob <- CellCycleScoring(ob, s.features = s.genes_ensembl, g2m.features = g2m.genes_ensembl, set.ident = TRUE)

ob <- ScaleData(ob, vars.to.regress = c("S.Score", "G2M.Score", "nUMI", "mitoRatio"), features = rownames(ob))



ob <- RunPCA(ob,features = VariableFeatures(object = ob))

#method in hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
# The first metric is the point where the PCs only contribute 5% of SD and the PCs cumulatively contribute 90% of the SD
#second metric ins the PC where the percent change in variation btwn consequtive PCs is less than 0.1%
#take minimum of the two metrics

pct = ob[["pca"]]@stdev / sum(ob[["pca"]]@stdev) * 100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct <5)[1]
co2 = sort(which((pct[1:length(pct) -1]-pct[2:length(pct)])>0.1),decreasing=T)[1] + 1
min(co1,co2)


##
ob <- FindNeighbors(ob, dims = 1:13)
ob <- FindClusters(ob, resolution = 0.5)
##

ob <- RunUMAP(ob, dims = 1:13)

DimPlot(ob, reduction = "umap", label=T, pt.size=0.86)



#can change only.pos to false if you want negative markers as well. for now i only want to see po\sitively diff. expressed markers in a cluster
ob.markers <- FindAllMarkers(ob, only.pos = TRUE)

save(ob, file = "./results/Rdata/seurat_ob.Rdata")
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters

FeaturePlot(ob, features = c("Rasd1"))
#markers distinguishing cluster

FeaturePlot(ob, features = c("Rasd1"), sort.cell = T, pt.size=0.86)



ob.markers[which(ob.markers$gene == "Qsox1"),] #cluster 7

#markers distinguishing cluster
cluster7.markers <- FindMarkers(ob, ident.1 = 7, min.pct = 0.25)




#########################################################################################




# 
# 
# 
# 
# #percent MT
# ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^mt-")
# 
# # Visualize QC metrics as a violin plot
# VlnPlot(ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# 
# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# 
# plot1 <- FeatureScatter(ob, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(ob, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))
# 
# ob <- subset(ob, subset = nFeature_RNA > 200 & nFeature_RNA < 5800 & percent.mt < 10)
# 
# 
# ob <- NormalizeData(ob, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# #identify highly variable features (high cell-to-cell variation)
# ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = 2000)
# 
# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(ob), 10)
# 
# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(ob)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# CombinePlots(plots = list(plot1, plot2))
# 
# #Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA
# all.genes <- rownames(ob)
# ob <- ScaleData(ob, features = all.genes)
# 
# ob <- RunPCA(ob,features = VariableFeatures(object = ob))
# 
# # Examine and visualize PCA results a few different ways
# print(ob[["pca"]], dims = 1:5, nfeatures = 5)
# 
# VizDimLoadings(ob, dims = 1:2, reduction = "pca")
# 
# DimPlot(ob, reduction = "pca")
# 
# #Determine the ‘dimensionality’ of the dataset
# 
# #jackstraw: randomly permute 1% of data, and rerun PCA, constructing a ‘null distribution’ of feature scores, and repeat this procedure. 
# #We identify ‘significant’ PCs as those who have a strong enrichment of low p-value features.
# 
# # ob <- JackStraw(ob, num.replicate = 100)
# # ob <- ScoreJackStraw(ob, dims = 1:20)
# # 
# # JackStrawPlot(ob, dims = 1:20) # useless
# 
# ElbowPlot(ob) #looks like crap
# 
# #try method in hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
# # The first metric is the point where the PCs only contribute 5% of SD and the PCs cumulatively contribute 90% of the SD
# #second metric ins the PC where the percent change in variation btwn consequtive PCs is less than 0.1%
# #take minimum of the two metrics
# pct = ob[["pca"]]@stdev / sum(ob[["pca"]]@stdev) * 100
# cumu = cumsum(pct)
# co1 = which(cumu > 90 & pct <5)[1]
# co2 = sort(which((pct[1:length(pct) -1]-pct[2:length(pct)])>0.1),decreasing=T)[1] + 1
# min(co1,co2)
# 
# 
# ##
# ob <- FindNeighbors(ob, dims = 1:13)
# ob <- FindClusters(ob, resolution = 0.5)
# ##
# 
# ob <- RunUMAP(ob, dims = 1:13)
# 
# DimPlot(ob, reduction = "umap", label=T)
# 
# # note that you can set `label = TRUE` or use the LabelClusters function to help label
# # individual clusters
# 
# FeaturePlot(ob, features = c("Rasd1"))
# 
# save(ob, file = "./results/Rdata/seurat_ob.Rdata")
# #load("./results/Rdata/seurat_ob.Rdata")
# #DimPlot(ob, reduction = "umap", label=TRUE,split.by = "seurat_clusters")
# 
# DimPlot(ob, reduction = "umap", label=T,pt.size = 0.86)
# 
# FeaturePlot(ob, features = c("Qsox1"), sort.cell = T, pt.size=0.86)
# 
# #cluster5.markers <- FindMarkers(ob, ident.1 = 0, min.pct = 0.25)
# 
# ob.markers <- FindAllMarkers(ob, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# ob.markers[which(ob.markers$gene == "Qsox1"),] #cluster 4
# #ob.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
# 
# #markers distinguishing cluster
# cluster4.markers <- FindMarkers(ob, ident.1 = 4, min.pct = 0.25)
# #cluster4.markers <- FindMarkers(ob, ident.1 = c(4), ident.2 = c(0,1,2,3,5,6,7,8,9,10), min.pct = 0.25)
# head(cluster4.markers, n = 30)
# 
# # FeaturePlot(ob, features = c("Qsox1"))
# # 
# # 
# # FeaturePlot(ob, features = c("Tnfsf11", "Tnfrsf11a", "Tnfrsf11b"))
# # 
# # 
# # FeaturePlot(ob, features = c("Sp7", "Ibsp", "Bglap", "Col1a1"))
# # 
# # FeaturePlot(ob, features = c("Sost", "Dmp1", "Phex", "Mepe"))
# # 
# # FeaturePlot(ob, features = c("Ctsk", "Ocstamp", "Dcstamp"))
# # FeaturePlot(ob, features = c("Adcy9"))
# # FeaturePlot(ob, features = c("Rasd1"))
