############
library(sva)
library(Mus.musculus)#seems to use mm10
library(WGCNA)
library(topGO)
library(DESeq2)
library("igraph")
library("bnlearn")
library("parallel")
library(rtracklayer)
############
options(stringsAsFactors = FALSE)
set.seed(8675309)
#

# read in the RNA-seq processed counts file
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)


#
annot_file = read.csv("~/Documents/projects/DO_project/results/flat/annot_file.csv")
annot_file = annot_file[,c(1,2)]

counts = counts[which(rownames(counts) %in% annot_file$Gene.ID),]
#

#find and remove features that have fewer than 10 reads in more than 90% (173) of samples 
x=c()
for(i in 1:nrow(counts)){
  if(sum(counts[i,]<10)>=173){
    print(i)
    x = append(x,i)
  }
}

#311 genes removed
counts = counts[-x,]



######


#vst from deseq2
vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))
c = as.numeric((unlist(strsplit(colnames(vst),"-"))))

colnames(vst) = c
#get batch. What I did here is I got batch from the file names of the alignment output for RNA-seq
#DOESNT TAKE RE-POOL INTO ACCOUNT, OR DIFFERENT RUN DATES FOR SPLIT POOLS
f = list.files("./results/flat/RNA-seq/sums/")
p1 = f[grep(pattern = "Pool1",f)]
p2 = f[grep(pattern = "Pool2",f)]
p3 = f[grep(pattern = "SetA",f)]
p4 = f[grep(pattern = "setB",f)]

b1 = c()
b2 = c()
b3 = c()
b4 = c()

for(i in 1:length(p1)){
  b1 = c(b1,strsplit(strsplit(p1[i],"Pool1_")[[1]][2],"-")[[1]][1])
  b1 = unique(b1)
}

for(i in 1:length(p2)){
  b2 = c(b2,strsplit(strsplit(p2[i],"Pool2_")[[1]][2],"-")[[1]][1])
  b2 = unique(b2)
}

for(i in 1:length(p3)){
  b3 = c(b3,strsplit(strsplit(strsplit(p3[i],"SetA_")[[1]][2],"_")[[1]][2],"-")[[1]][1])
  b3 = unique(b3)
}

for(i in 1:length(p4)){
  b4 = c(b4,strsplit(strsplit(strsplit(p4[i],"setB_")[[1]][2],"_")[[1]][2],"-")[[1]][1])
  b4 = unique(b4)
}


#read in the raw phenotypes file
x = read.csv("./results/flat/full_pheno_table.csv", stringsAsFactors = FALSE)

#get covs by matching colnames of vst with mouse ID in raw pheno file
covs = x[match(colnames(vst),x$Mouse.ID),]

#sac.date,sex, age at sac, generation
covs = covs[,c(2,5,7,22)]

covs$male = NA
covs$female = NA
covs$batch1=NA
covs$batch2=NA
covs$batch3=NA
covs$batch4=NA



#define batches
covs$batch1[match(b1,rownames(covs))] = 1
covs[which(rownames(covs)%in% b1==FALSE),"batch1"] = 0

covs$batch2[match(b2,rownames(covs))] = 1
covs[which(rownames(covs)%in% b2==FALSE),"batch2"] = 0

covs$batch3[match(b3,rownames(covs))] = 1
covs[which(rownames(covs)%in% b3==FALSE),"batch3"] = 0

covs$batch4[match(b4,rownames(covs))] = 1
covs[which(rownames(covs)%in% b4==FALSE),"batch4"] = 0



covs = as.data.frame(covs[,2:10])

covs$batch = NA

covs[which(covs$batch1 == 1),"batch"] = 1
covs[which(covs$batch2 == 1),"batch"] = 2
covs[which(covs$batch3 == 1),"batch"] = 3
covs[which(covs$batch4 == 1),"batch"] = 4


batch = covs$batch

#generation is confounded with batch so it is not added
modcombat = model.matrix(~as.factor(sex) + as.factor(age_at_sac_days), data=covs)

#batch removal
edata = ComBat(dat=vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#transpose the matrix
edata = t(edata)

save(edata, file = "./results/Rdata/networks/edata_full.Rdata")

#check data
gsg = goodSamplesGenes(edata, verbose = 3);
gsg$allOK

##
#cluster samples
sampleTree = hclust(dist(edata), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#####



#pick the soft thresholding power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(edata, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE)
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
###


#make network with sft=4
net = blockwiseModules(edata, power = 4,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE,
                       saveTOMs = FALSE,
                       verbose = 3)

## 4: 39 modules not including 0, 6004 genes in module 0

#saveRDS(net, file="./results/Rdata/networks/wgcna_4.RDS")
net = readRDS("./results/Rdata/networks/wgcna_4.RDS")

#get traits we want to look at
pheno = read.csv("./results/flat/full_pheno_table.csv", stringsAsFactors = FALSE)


datTraits = pheno[which(pheno$Mouse.ID %in% rownames(edata)),]
datTraits = datTraits[match(rownames(edata),datTraits$Mouse.ID),]

#remove non-pheno columns
datTraits = datTraits[,-c(1:7,16,22)]

#convert to numeric
for(i in 1:ncol(datTraits)){datTraits[,i] = as.numeric(datTraits[,i])}

#remove histo columns 
#datTraits = datTraits[,-c(26:43)]
#remove adipose
#datTraits = datTraits[,-c(44:55)]

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)


# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Construct numerical labels corresponding to the colors
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
nGenes = ncol(edata);
nSamples = nrow(edata);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(edata, moduleColors)$eigengenes
MEs = orderMEs(MEs0)



#REMOVE GREY
MEs = MEs[,-which(colnames(MEs) == "MEgrey")]
#cor module eigengenes with traits
moduleTraitCor = cor(MEs, datTraits, use = "p",method = "s");

moduleTraitPvalue = as.data.frame(matrix(nrow = nrow(moduleTraitCor),ncol = ncol(moduleTraitCor)))
for(i in 1:ncol(moduleTraitCor)){
  nSamples = length(which(is.na(datTraits[,colnames(moduleTraitCor)[i]]) == FALSE))
  moduleTraitPvalue[,i] = corPvalueStudent(moduleTraitCor[,i], nSamples) # uses sample for each trait. is this correct?
  print(colnames(moduleTraitCor)[i])
  print(nSamples)
  
}
nSamples=nrow(edata)
colnames(moduleTraitPvalue) = colnames(moduleTraitCor)
rownames(moduleTraitPvalue) = rownames(moduleTraitCor)

##remove everything but uCT
#moduleTraitPvalue = moduleTraitPvalue[,c(44:61)]
#moduleTraitCor = moduleTraitCor[,c(44:61)]

#moduleTraitPvalue = moduleTraitPvalue[,c(11:61)]
#moduleTraitCor = moduleTraitCor[,c(11:61)]

#moduleTraitPvalue = as.matrix(moduleTraitPvalue)

#keep all but  weight, length, glucose , fat pads, muscle masses and MAT
moduleTraitPvalue = moduleTraitPvalue[,c(11:61)]
moduleTraitCor = moduleTraitCor[,c(11:61)]


moduleTraitPvalue = as.matrix(moduleTraitPvalue)

save(moduleTraitPvalue, file = "./results/Rdata/networks/moduleTraitPvalue_full_4.RData")
save(moduleTraitCor, file = "./results/Rdata/networks/moduleTraitCor_full_4.RData")

#identify significant modules
sig_mod = moduleTraitPvalue[which(rownames(moduleTraitPvalue) %in% names(which(apply(moduleTraitPvalue, 1, function(r) any(r < 0.05/ncol(MEs)))))),]

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = signif(moduleTraitPvalue, 1)

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

trait_names = gsub(pattern = "bending_",replacement = "",x = names(datTraits))
trait_names = gsub(pattern = "uCT_",replacement = "",x = trait_names)
trait_names = gsub(pattern = "\\.\\.",replacement = "\\.",x = trait_names)
trait_names = gsub(pattern = "\\.\\.",replacement = "\\.",x = trait_names)
trait_names[8] = "Adiposity"

####
####
#only trabecular traits
moduleTraitCor = moduleTraitCor[,c(34:41)]

trait_names = colnames(moduleTraitCor)
trait_names = gsub(pattern = "uCT_",replacement = "",x = trait_names)
#MEs = MEs[,-37]
textMatrix = signif(moduleTraitPvalue, 1)
####
####

#dim(textMatrix) = dim(moduleTraitCor)
textMatrix = textMatrix[,c(34:41)]
# Display the correlation values within a heatmap plot
par(mar=c(4.1, 13.1, 3.1, 2.1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = trait_names,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.6,
               zlim = c(-1,1),
               cex.lab.x = 0.85,
               cex.lab.y = 0.7,
               verticalSeparator.x = c(1:length(trait_names)),
               horizontalSeparator.y = c(1:length(names(MEs))),
               main = paste("Module-trait relationships"))
###
which(moduleTraitCor == sort((moduleTraitCor),decreasing =TRUE)[1], arr.ind = T)
colnames(moduleTraitCor)[33]


#create trimmed expression object to remove grey module and get gene names from annot file
########################
combat_annot = as.data.frame(colnames(edata))

combat_annot$module = net$colors
combat_annot$color = moduleColors

#REMOVE GREY
rmv = combat_annot[which(combat_annot$color == "grey"),"colnames(edata)"]

combat_annot = combat_annot[-which(combat_annot$`colnames(edata)` %in% rmv),]
edata_trim = edata[,-(which(colnames(edata) %in% rmv))]

combat_annot = combat_annot[,-2]

combat_annot[,c(3:4)] = annot_file[match(combat_annot$`colnames(edata)`,annot_file$Gene.ID),c(1,2)]


geneModuleMembership = as.data.frame(cor(edata_trim, MEs, use = "p")) # spearman or kendall? using pearson here

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership$gene = colnames(edata_trim)


combat_annot[5:(ncol(geneModuleMembership)+4)] = geneModuleMembership[match(geneModuleMembership$gene,combat_annot$`colnames(edata)`),]


save(combat_annot, file = "./results/Rdata/networks/geneModMemAnnot_power4.RData")

save(moduleTraitPvalue, file = "./results/Rdata/networks/moduleTraitPvalue_full_4.RData")

save(moduleTraitCor, file = "./results/Rdata/networks/moduleTraitCor_full_4.RData")



#Do GO analysis
##
load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
combat_annot = geneModMemAnnot
moduleColors = moduleColors[-which(moduleColors=="grey")]
hubs = chooseTopHubInEachModule(edata,moduleColors)


network_GO = list()
for(color in unique(combat_annot$color)){
  
  allGenes = combat_annot$gene
  interesting.genes = combat_annot[which(combat_annot$color == color),"gene"]

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
  network_GO[[color]] = t.all
}
save(network_GO, file="./results/Rdata/networks/GO_sft4.RData")
#####

