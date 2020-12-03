
library(clustree)
library(Seurat)
library(topGO)
library(Hmisc)
set.seed(8675309)


#Read in the data
ob.data <- Read10X(data.dir = "./data/filtered_gene_bc_matrices/mm10/",unique.features = T)




#get cell cycle markers (S and G2/M, from Tirosh et al, 2015) genes, put them in the format Aa*
  

s.genes <- cc.genes$s.genes

s.genes<-tolower(s.genes)

s.genes<-capitalize(s.genes)


g2m.genes <- cc.genes$g2m.genes

g2m.genes<-tolower(g2m.genes)

g2m.genes<-capitalize(g2m.genes)



#create Seurat object
#only use features (genes) detected in at leasr 3 cells, and include cells where at least 200 features are detected. Filters out a few empty droplets probably


ob <- CreateSeuratObject(counts = ob.data, project = "scRNA", min.cells = 3, min.features = 200)


#PREPROCESSING
#find percentage of counts in a cell coming from mitochondrial genes. High mito contamination usually means low quality/dying cells.

ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^mt-")

#Plot QC parameters
VlnPlot(ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#Based on plot, only take with more than 800 reads (removes that bottom one with 719) and less than 5800 (removes those three outliers at the top). 
#Also remove anything with more than 10% mito reads


ob <- subset(ob, subset = nFeature_RNA > 800 & nFeature_RNA < 5800 & percent.mt < 10)



#Normalize feature expression measurements for each cell by total expression, multiply by 10,000 and log normalize, then identify 3000 most highly variable features in dataset, focus on them in downstream analysis.
#Finally, scale the data (linear transformation)



ob <- NormalizeData(ob)
ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = 3000)
ob <- ScaleData(ob, features = rownames(ob))


#Identify top 10 most variable and plot them


top10 = head(VariableFeatures(ob),10)
plot1 = VariableFeaturePlot(ob)
plot2 = LabelPoints(plot=plot1, points = top10, repel=T)             

plot2



#identify sex

#Xist_expression = (FetchData(ob, "Xist"))

#if cell has 0% xist expression, make male. else, female
ob[["percent.xist"]] <- PercentageFeatureSet(ob, pattern = "Xist")

ob[["sex"]] = 0

ob[["sex"]][ob[["percent.xist"]] == 0] = 1


#Run PCA using variable genes


ob = RunPCA(ob, features = VariableFeatures(ob))


#score cells by cell cycle genes

ob <- CellCycleScoring(ob, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#Run PCA on the cell cycle genes and plot

ob <- RunPCA(ob, features = c(s.genes, g2m.genes), approx=FALSE)
DimPlot(ob)

#Cells cluster by phase, but theres also something else in the data with a huge effect in PC1 and also PC2.

#Regress out cell cycle scores and percentMT too, then run PCA on the variable genes ID'd previously and plot

ob <- ScaleData(ob, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(ob))


# Run PCA on variable features after regression
ob <- RunPCA(ob, features = VariableFeatures(ob), nfeatures.print = 10)
DimPlot(ob)



#still some weird variation in PC2. Is it sex?
#Plot and color by sex
DimPlot(ob, group.by = "sex")


#Its not sex. immune-related/myeloid cells/osteoclasts? Look at top negative PC2 genes. To me, top negative PC2's appear to be immune related.
#Here's some plots of the top negative PC2 genes


par(mfrow=c(2,2))
FeaturePlot(ob, features = "Tyrobp")
FeaturePlot(ob, features = "Lyz2")
FeaturePlot(ob, features = "Trem2")
FeaturePlot(ob, features = "Ctss")




#Code here determines ow many PCs to use. From hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html

# Determine percent of variation associated with each PC
pct <- ob[["pca"]]@stdev / sum(ob[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
#13


#Use 13 PCs.


##clustering 
ob <- FindNeighbors(object = ob, dims=1:13)



#What resolution to use in clustering? We can use clustree

###### choose resolution
ob_t <- FindClusters(object = ob, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5))

clustree(ob_t, prefix= "RNA_snn_res.", node_colour = "sc3_stability")

#I chose 1.0. Find clusters and run UMAP


ob <- FindClusters(object = ob, resolution = 1)

ob <- RunUMAP(ob, dims = 1:13)




DimPlot(ob,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        plot.title = "UMAP")




##remove the far away cluster, UMAP1 close to 15. 13 total cells
outliers = names(which(ob@reductions$umap@cell.embeddings[,"UMAP_1"] > 14))
not_outliers = names(which(ob@reductions$umap@cell.embeddings[,"UMAP_1"] < 14))
ob_s = subset(ob, cells = not_outliers)
#ID all cluster markers, also clusters for genes we're interested in.



DimPlot(ob_s,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        plot.title = "UMAP")

ob.markers <- FindAllMarkers(ob, only.pos = TRUE)


save(ob, file = "./results/Rdata/seurat_ob.Rdata")


# #First, Rasd1 markers. (Cluster 10)
# 
# ob.markers[which(ob.markers$gene =="Rasd1"),] # 10
# FeaturePlot(ob, features = "Rasd1", sort.cell = T, pt.size = 1)
# 
# #output the top 30 cluster 10 genes
# 
# cluster10.markers <- FindMarkers(ob, ident.1 = 10, min.pct = 0.25,)
# 
# head(cluster10.markers, 30)




#Then, Qsox1 markers. (Cluster 1)

ob.markers[which(ob.markers$gene =="Qsox1"),] #cluster 1

#FeaturePlot(ob, features = "Qsox1", sort.cell = T, pt.size = 1)


cluster1.markers <- FindMarkers(ob, ident.1 = 1, min.pct = 0.25)

#Top cluster 1 genes

head(cluster1.markers, 30)





##GO analysis, but only for positive genes in a cluster
 

allGenes = rownames(ob)

cluster1.markers.pos <- FindMarkers(ob, ident.1 = 1, min.pct = 0.25,only.pos = T)
cluster10.markers.pos <- FindMarkers(ob, ident.1 = 10, min.pct = 0.25, only.pos = T)

interesting.genes = rownames(cluster1.markers.pos)

geneList<-factor(as.integer(allGenes %in% interesting.genes)) #If TRUE returns 1 as factor, otherwise 0
names(geneList)<-allGenes
###MF###
GOdata <- new("topGOdata", ontology = "MF", allGenes =geneList,
              annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t1<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t1)
###CC###
GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
              annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t2<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t2)
###BP###
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t3<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t3)
####
t.all = NULL
t.all<-rbind(t1,t2,t3)
t.all$classic<-as.numeric(as.character(t.all$classic))
#####









