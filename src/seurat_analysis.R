library(clustree)
library(Seurat)
library(topGO)
library(Hmisc)
set.seed(8675309)

#read in the data
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


##PREPROCESSING

#find percentage of counts in a cell coming from mitochondrial genes. High mito contamination usually means low quality/dying cells.
ob[["percent.mt"]] <- PercentageFeatureSet(ob, pattern = "^mt-")

#plot
VlnPlot(ob, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#Based on plot, only take with more than 800 reads (removes that bottom one with 719) and less than 5800 (removes those three outliers at the top). 
#Removes anythin with more than 10% mito reads
ob <- subset(ob, subset = nFeature_RNA > 800 & nFeature_RNA < 5800 & percent.mt < 10)



#normalize feature expression measurements for each cell by total expression, multiply by 10,000 and log normalize
ob <- NormalizeData(ob)

#Identify highly variable features in dataset, focus on them in downstream analysis. 2000 features returned by default
ob <- FindVariableFeatures(ob, selection.method = "vst")

#identify top 10 more variable and plot them
top10 = head(VariableFeatures(ob),10)
plot1 = VariableFeaturePlot(ob)
plot2 = LabelPoints(plot=plot1, points = top10, repel=T)             

plot2


#scale the data. linear transformation
ob <- ScaleData(ob, features = rownames(ob))



#identify sex

Xist_epression = (FetchData(ob, "Xist"))

#if cell has 0% xist expression, make male. else, female
ob[["percent.xist"]] <- PercentageFeatureSet(ob, pattern = "Xist")

ob[["sex"]] = 0

ob[["sex"]][ob[["percent.xist"]] == 0] = 1

#run PCA using variable genes

ob = RunPCA(ob, features = VariableFeatures(ob))


#score cells by cell cycle genes
ob <- CellCycleScoring(ob, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


#RunPCA on the cell cycle genes
ob <- RunPCA(ob, features = c(s.genes, g2m.genes), approx=FALSE)
DimPlot(ob)
#cells cluster by phase, but theres also something else in the data with a huge effect in PC1 and also PC2. Might be sex?


#regress out cell cycle scores and percentMT too
ob <- ScaleData(ob, vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"), features = rownames(ob))


# Run PCA on variable features after regression
ob <- RunPCA(ob, features = VariableFeatures(ob), nfeatures.print = 10)
DimPlot(ob)



# #still some weirf variation in PC2. try regressing out sex as well
# ob_s <- ScaleData(ob, vars.to.regress = c("S.Score","G2M.Score", "percent.mt","sex"), features = rownames(ob))
# 
# ob_s <- RunPCA(ob_s, features = VariableFeatures(ob_s), nfeatures.print = 10)
# DimPlot(ob_s, group.by = "sex")

##its not sex. immune-related/myeloid cells/osteoclasts?
#to me, top negative PC2's appear to be immune related
FeaturePlot(ob, features = "Tyrobp")




##determine PCA params
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
#12



##clustering 
ob <- FindNeighbors(object = ob, dims=1:12)

###### choose resolution
ob_t <- FindClusters(object = ob, resolution = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5))#0.8?

clustree(ob_t, prefix= "RNA_snn_res.", node_colour = "sc3_stability")
######
ob <- FindClusters(object = ob, resolution = 0.8)

ob <- RunUMAP(ob, dims = 1:12)





DimPlot(ob,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        plot.title = "UMAP")
dev.off()#



#ID al markers, also clusters for genes we're interested in
ob.markers <- FindAllMarkers(ob, only.pos = TRUE)

#save(ob, file = "./results/Rdata/seurat_ob.Rdata")


ob.markers[which(ob.markers$gene =="Rasd1"),] # 8
FeaturePlot(ob, features = "Rasd1", sort.cell = T, pt.size = 1)
#markers distinguishing cluster
#
cluster8.markers <- FindMarkers(ob, ident.1 = 8, min.pct = 0.25,)


ob.markers[which(ob.markers$gene =="Qsox1"),] #cluster 5

FeaturePlot(ob, features = "Qsox1", sort.cell = T, pt.size = 1)

cluster5.markers <- FindMarkers(ob, ident.1 = 5, min.pct = 0.25)



##GO analysis, but only for positive genes in a cluster
allGenes = rownames(ob)

cluster5.markers.pos <- FindMarkers(ob, ident.1 = 5, min.pct = 0.25,only.pos = T)
cluster8.markers.pos <- FindMarkers(ob, ident.1 = 8, min.pct = 0.25, only.pos = T)

interesting.genes = rownames(cluster5.markers.pos)

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
######
######
######









