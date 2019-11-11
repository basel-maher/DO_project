##wgcna sex specific  
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
#

#read in VST transformed, quantile-normalized data (produced in ./src/normalize_RNAseq.R)
#load("./results/Rdata/counts_vst_qnorm.Rdata")
#For now try non-quantile normalized data, compare with q-normed

# read in the RNA-seq processed counts file
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#read in an annotation file. This is output from the RNA-seq pipeline
#annot_file = read.delim("./data/314-FarberDO2_S6.gene_abund.tab",header = TRUE) #this looks crappy
annot_file = readGFF("./results/flat/RNA-seq/mus_stringtie_merged.gtf")
annot_file = annot_file[which(annot_file$type=="transcript"),]
#remove transcripts not on somatic and X chroms
chr = c(seq(1:19),"X")
annot_file = annot_file[which(annot_file$seqid %in% chr),]


annot_file = annot_file[,c(9,11)]
annot_file = unique(annot_file)
#

counts = counts[which(rownames(counts) %in% annot_file$gene_id),]
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


#some genes are duplicated
counts = t(counts)

colnames(counts) = annot_file[match(colnames(counts),annot_file$gene_id),"gene_name"]


#for now, give them a unique ID
colnames(counts)[which(duplicated(colnames(counts)))] = paste0(colnames(counts)[which(duplicated(colnames(counts)))], "_isoform")

which(duplicated(colnames(counts)))
colnames(counts)[3570] = paste0(colnames(counts)[3570],".2")
colnames(counts)[5342] = paste0(colnames(counts)[5342],".2")
colnames(counts)[6661] = paste0(colnames(counts)[6661],".3")

#colnames(counts)[4133] = paste0(colnames(counts)[4133],".2")
#colnames(counts)[4134] = paste0(colnames(counts)[4134],".3")
#colnames(counts)[18694] = paste0(colnames(counts)[18694],".2")

counts = t(counts)


######


#vst from deseq2
vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))
c = as.numeric((unlist(strsplit(colnames(vst),"-"))))
#c = c[-which(is.na(c))]
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


covs$batch1=NA
covs$batch2=NA
covs$batch3=NA
covs$batch4=NA



#covs$male = as.numeric(covs$sex=="M")
#covs$female = as.numeric(covs$sex=="F")
#define batches
covs$batch1[match(b1,rownames(covs))] = 1
covs[which(rownames(covs)%in% b1==FALSE),"batch1"] = 0

covs$batch2[match(b2,rownames(covs))] = 1
covs[which(rownames(covs)%in% b2==FALSE),"batch2"] = 0

covs$batch3[match(b3,rownames(covs))] = 1
covs[which(rownames(covs)%in% b3==FALSE),"batch3"] = 0

covs$batch4[match(b4,rownames(covs))] = 1
covs[which(rownames(covs)%in% b4==FALSE),"batch4"] = 0


#cc = as.matrix(covs[,c(3,5:8)])
#rownames(cc) = rownames(covs)

#convert to integers
#covs[,1] = as.factor(covs[,1])

#cc[,c(1:5)] = as.integer(cc[,])

covs = as.data.frame(covs[,2:10])

covs$batch = NA

covs[which(covs$batch1 == 1),"batch"] = 1
covs[which(covs$batch2 == 1),"batch"] = 2
covs[which(covs$batch3 == 1),"batch"] = 3
covs[which(covs$batch4 == 1),"batch"] = 4


batch = covs$batch

#generation is confounded with batch so it is not added
modcombat = model.matrix(~as.factor(age_at_sac_days), data=covs)

#batch removal
edata = ComBat(dat=vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#transpose the matrix
edata = t(edata)


#check data
gsg = goodSamplesGenes(edata, verbose = 3);
gsg$allOK

males = x[which(x$sex=="M"),"Mouse.ID"]
edata_m = edata[which(rownames(edata) %in% males),]
edata_f = edata[which(rownames(edata) %in% males == FALSE),]
##


#cluster samples
sampleTree_m = hclust(dist(edata_m), method = "average")
sampleTree_f = hclust(dist(edata_f), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree_f, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#####



#pick the soft thresholding power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft_m = pickSoftThreshold(edata_m, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE)
sft_f = pickSoftThreshold(edata_f, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft_m$fitIndices[,1], -sign(sft_m$fitIndices[,3])*sft_m$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft_m$fitIndices[,1], -sign(sft_m$fitIndices[,3])*sft_m$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft_m$fitIndices[,1], sft_m$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_m$fitIndices[,1], sft_m$fitIndices[,5], labels=powers, cex=cex1,col="red")
###
#do 5 for males
net_m = blockwiseModules(edata_m, power = 5,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)

##7: 27 mods no 0, 4888 in mod  0

#try 4 for females, 4 is the 0.9 thresh
net_f = blockwiseModules(edata_f, power = 4,
                         TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = FALSE,
                         verbose = 3)

##7: 17 mods no 0, 8357 in in mod 0
#
#
#
#
#get traits we want to look at
pheno = read.csv("./results/flat/full_pheno_table.csv", stringsAsFactors = FALSE)


datTraits_m = pheno[which(pheno$Mouse.ID %in% rownames(edata_m)),]
datTraits_f = pheno[which(pheno$Mouse.ID %in% rownames(edata_f)),]

datTraits_m = datTraits_m[match(rownames(edata_m),datTraits_m$Mouse.ID),]
datTraits_f = datTraits_f[match(rownames(edata_f),datTraits_f$Mouse.ID),]

#remove non-pheno columns
datTraits_m = datTraits_m[,-c(1:7,16,22)]
datTraits_f = datTraits_f[,-c(1:7,16,22)]

#convert to numeric
for(i in 1:ncol(datTraits_m)){datTraits_m[,i] = as.numeric(datTraits_m[,i])}
for(i in 1:ncol(datTraits_f)){datTraits_f[,i] = as.numeric(datTraits_f[,i])}

#remove histo columns 
#datTraits = datTraits[,-c(26:43)]
#remove adipose
#datTraits = datTraits[,-c(44:55)]

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors_f = labels2colors(net_f$colors)
mergedColors_m = labels2colors(net_m$colors)


# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_m$dendrograms[[1]], mergedColors[net_m$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Construct numerical labels corresponding to the colors
moduleLabels_m = net_m$colors
moduleLabels_f = net_f$colors

moduleColors_m = labels2colors(net_m$colors)
moduleColors_f = labels2colors(net_f$colors)

MEs_m = net_m$MEs;
MEs_f = net_f$MEs;

geneTree_m = net_m$dendrograms[[1]]
geneTree_f = net_f$dendrograms[[1]]

nGenes_m = ncol(edata_m);
nGenes_f = ncol(edata_f);

nSamples_m = nrow(edata_m);
nSamples_f = nrow(edata_f);

# Recalculate MEs with color labels
MEs0_m = moduleEigengenes(edata_m, moduleColors_m)$eigengenes
MEs0_f = moduleEigengenes(edata_f, moduleColors_f)$eigengenes

MEs_m = orderMEs(MEs0_m)
MEs_f = orderMEs(MEs0_f)


#REMOVE GREY
MEs_m = MEs_m[,-which(colnames(MEs_m) == "MEgrey")]
MEs_f = MEs_f[,-which(colnames(MEs_f) == "MEgrey")]

#cor module eigengenes with traits
moduleTraitCor_m = cor(MEs_m, datTraits_m, use = "p",method = "s");
moduleTraitCor_f = cor(MEs_f, datTraits_f, use = "p",method = "s");

moduleTraitPvalue_m = as.data.frame(matrix(nrow = nrow(moduleTraitCor_m),ncol = ncol(moduleTraitCor_m)))

for(i in 1:ncol(moduleTraitCor_m)){
  nSamples_m = length(which(is.na(datTraits_m[,colnames(moduleTraitCor_m)[i]]) == FALSE))
  moduleTraitPvalue_m[,i] = corPvalueStudent(moduleTraitCor_m[,i], nSamples_m) # uses sample for each trait. is this correct?
  print(colnames(moduleTraitCor_m)[i])
  print(nSamples_m)
  
}
nSamples_m=nrow(edata_m)
colnames(moduleTraitPvalue_m) = colnames(moduleTraitCor_m)
rownames(moduleTraitPvalue_m) = rownames(moduleTraitCor_m)
#
#
#
moduleTraitPvalue_f = as.data.frame(matrix(nrow = nrow(moduleTraitCor_f),ncol = ncol(moduleTraitCor_f)))

for(i in 1:ncol(moduleTraitCor_f)){
  nSamples_f = length(which(is.na(datTraits_f[,colnames(moduleTraitCor_f)[i]]) == FALSE))
  moduleTraitPvalue_f[,i] = corPvalueStudent(moduleTraitCor_f[,i], nSamples_f) # uses sample for each trait. is this correct?
  print(colnames(moduleTraitCor_f)[i])
  print(nSamples_f)
  
}
nSamples_f=nrow(edata_f)
colnames(moduleTraitPvalue_f) = colnames(moduleTraitCor_f)
rownames(moduleTraitPvalue_f) = rownames(moduleTraitCor_f)

##remove everything but uCT
#moduleTraitPvalue = moduleTraitPvalue[,c(44:61)]
#moduleTraitCor = moduleTraitCor[,c(44:61)]

#moduleTraitPvalue = moduleTraitPvalue[,c(11:61)]
#moduleTraitCor = moduleTraitCor[,c(11:61)]

#moduleTraitPvalue = as.matrix(moduleTraitPvalue)

#keep all but  weight, length, glucose , fat pads, muscle masses and MAT
moduleTraitPvalue_m = moduleTraitPvalue_m[,c(11:61)]
moduleTraitCor_m = moduleTraitCor_m[,c(11:61)]

moduleTraitPvalue_f = moduleTraitPvalue_f[,c(11:61)]
moduleTraitCor_f = moduleTraitCor_f[,c(11:61)]


moduleTraitPvalue_m = as.matrix(moduleTraitPvalue_m)
moduleTraitPvalue_f = as.matrix(moduleTraitPvalue_f)


sig_mod_m = moduleTraitPvalue_m[which(rownames(moduleTraitPvalue_m) %in% names(which(apply(moduleTraitPvalue_m, 1, function(r) any(r < 0.05/ncol(MEs_m)))))),]
sig_mod_f = moduleTraitPvalue_f[which(rownames(moduleTraitPvalue_f) %in% names(which(apply(moduleTraitPvalue_f, 1, function(r) any(r < 0.05/ncol(MEs_f)))))),]

######
#####
#####
#####
######
combat_annot_m = as.data.frame(colnames(edata_m))

combat_annot_m$module = net_m$colors
combat_annot_m$color = moduleColors_m

#REMOVE GREY
rmv = combat_annot_m[which(combat_annot_m$color == "grey"),"colnames(edata_m)"]

combat_annot_m = combat_annot_m[-which(combat_annot_m$`colnames(edata_m)` %in% rmv),]
edata_trim_m = edata_m[,-(which(colnames(edata_m) %in% rmv))]

#have to run this sometimes?
combat_annot_m = combat_annot_m[,-2]
#the gsub allows for matching of genes that had _isoform* added to them
combat_annot_m[,c(3:4)] = annot_file[match(gsub(combat_annot_m$`colnames(edata_m)`,pattern = "_isoform.*",replacement = ""),annot_file$gene_name),c(1,2)]


#modNames = substring(names(MEs), 3)
geneModuleMembership_m = as.data.frame(cor(edata_trim_m, MEs_m, use = "p")) # spearman or kendall? using pearson here

MMPvalue_m = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_m), nSamples_m))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership_m$gene = colnames(edata_trim_m)
#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="");


combat_annot_m[5:(ncol(geneModuleMembership_m)+4)] = geneModuleMembership_m[match(geneModuleMembership_m$gene,combat_annot_m$`colnames(edata_m)`),]
#
#
#
#
combat_annot_f = as.data.frame(colnames(edata_f))

combat_annot_f$module = net_f$colors
combat_annot_f$color = moduleColors_f

#REMOVE GREY
rmv = combat_annot_f[which(combat_annot_f$color == "grey"),"colnames(edata_f)"]

combat_annot_f = combat_annot_f[-which(combat_annot_f$`colnames(edata_f)` %in% rmv),]
edata_trim_f = edata_f[,-(which(colnames(edata_f) %in% rmv))]

#have to run this sometimes?
combat_annot_f = combat_annot_f[,-2]
#the gsub allows for matching of genes that had _isoform* added to them
combat_annot_f[,c(3:4)] = annot_file[match(gsub(combat_annot_f$`colnames(edata_f)`,pattern = "_isoform.*",replacement = ""),annot_file$gene_name),c(1,2)]


#modNames = substring(names(MEs), 3)
geneModuleMembership_f = as.data.frame(cor(edata_trim_f, MEs_f, use = "p")) # spearman or kendall? using pearson here

MMPvalue_f = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_f), nSamples_f))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership_f$gene = colnames(edata_trim_f)
#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="");


combat_annot_f[5:(ncol(geneModuleMembership_f)+4)] = geneModuleMembership_f[match(geneModuleMembership_f$gene,combat_annot_f$`colnames(edata_f)`),]







#########
#########
allGenes = combat_annot_f$gene
interesting.genes = combat_annot_f[which(combat_annot_f$color == "turquoise"),"gene"]
#blue is bone module, 2386 genes
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
t.all<-subset(t.all,t.all$classic<=0.01)
t.all<-t.all[order(t.all$classic,decreasing=FALSE),]
dim(t.all[t.all$classic<=1e-5,])
######

##############CONSENSUS###############
#put both expression datasets in one list
nSets = 2
setLabels = c("females", "males")
multiExpr = vector(mode = "list", length = nSets)
#rows = samples, cols = genes
multiExpr[[1]] = list(data = as.data.frame(edata_f))
multiExpr[[2]] = list(data = as.data.frame(edata_m))
#check
exprSize = checkSets(multiExpr)
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK

#cluster samples separately in each set
sampleTrees = list()
for (set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),xlab="", sub="", cex = 0.7)
}
####


#no outliers (follow outlier removal based on but height if u see anything. Horvath tutorial page)
#PHENO DATA
pheno = read.csv("./results/flat/full_pheno_table.csv", stringsAsFactors = FALSE)

datTraits = pheno[,-c(2:7,16,22)]

Traits = vector(mode="list", length = nSets);

for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, datTraits$Mouse.ID);
  Traits[[set]] = list(data = datTraits[traitRows, -1]);
  rownames(Traits[[set]]$data) = datTraits[traitRows, 1];
  }
collectGarbage()


# Define data set dimensions
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples


#######
# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}

#try 8

softPower = 8;
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate adjacencies in each individual data set
for (set in 1:nSets)adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower

# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes));
# Calculate TOMs in each individual data set
for (set in 1:nSets)TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])

##DONE ON RIVANNA HPC
###scale TOMs
#scale so that 95th percentiles of each sex match
scaleP=0.95
set.seed(8675309)
#sample TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000);
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list()
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)

# Loop over sets
for (set in 1:nSets)
{
  # Select the sampled TOM entries
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  # Calculate the 95th percentile
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  # Scale the male TOM
  if (set>1)
  {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}
#####
load("~/Desktop/TOM.Rdata")


#calculate consensus TOM
consensusTOM = pmin(TOM[1, , ], TOM[2, , ])

# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,deepSplit = 2, cutHeight = 0.995,minClusterSize = minModuleSize,pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)

save(TOM, file="./TOM.Rdata")