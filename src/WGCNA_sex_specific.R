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
set.seed(8675309)

# read in the RNA-seq processed counts file
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

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



covs = as.data.frame(covs[,2:8])

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
save(edata_m, file = "./results/Rdata/networks/edata_m.Rdata")
edata_f = edata[which(rownames(edata) %in% males == FALSE),]
save(edata_f, file = "./results/Rdata/networks/edata_f.Rdata")




#cluster samples
sampleTree_m = hclust(dist(edata_m), method = "average")
sampleTree_f = hclust(dist(edata_f), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree_m, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
plot(sampleTree_f, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

#####



#pick the soft thresholding power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft_m = pickSoftThreshold(edata_m, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE)


sft_f = pickSoftThreshold(edata_f, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE)


# Plot the results (do for M and F):
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

#SFT = 4 for females
#SFT = 5 for males


net_m = blockwiseModules(edata_m, power = 5,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE,
                       saveTOMs = FALSE,
                       verbose = 3)


saveRDS(net_m, file="./results/Rdata/networks/wgcna_m_5.RDS")

net_f = blockwiseModules(edata_f, power = 4,
                         TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.15,
                         numericLabels = TRUE,
                         saveTOMs = FALSE,
                         verbose = 3)

saveRDS(net_f, file="./results/Rdata/networks/wgcna_f_4.RDS")



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
plotDendroAndColors(net_f$dendrograms[[1]], mergedColors_f[net_f$blockGenes[[1]]],
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


#keep all but  weight, length, glucose , fat pads, muscle masses 
moduleTraitPvalue_m = moduleTraitPvalue_m[,c(11:65)]
moduleTraitCor_m = moduleTraitCor_m[,c(11:65)]

moduleTraitPvalue_f = moduleTraitPvalue_f[,c(11:65)]
moduleTraitCor_f = moduleTraitCor_f[,c(11:65)]


moduleTraitPvalue_m = as.matrix(moduleTraitPvalue_m)
moduleTraitPvalue_f = as.matrix(moduleTraitPvalue_f)


sig_mod_m = moduleTraitPvalue_m[which(rownames(moduleTraitPvalue_m) %in% names(which(apply(moduleTraitPvalue_m, 1, function(r) any(r < 0.05/ncol(MEs_m)))))),]
sig_mod_f = moduleTraitPvalue_f[which(rownames(moduleTraitPvalue_f) %in% names(which(apply(moduleTraitPvalue_f, 1, function(r) any(r < 0.05/ncol(MEs_f)))))),]


#create trimmed expression object to remove grey module and get gene names from annot file
#
#
#
#

combat_annot_m = as.data.frame(colnames(edata_m))

combat_annot_m$module = net_m$colors
combat_annot_m$color = moduleColors_m

#REMOVE GREY
rmv = combat_annot_m[which(combat_annot_m$color == "grey"),"colnames(edata_m)"]

combat_annot_m = combat_annot_m[-which(combat_annot_m$`colnames(edata_m)` %in% rmv),]
edata_trim_m = edata_m[,-(which(colnames(edata_m) %in% rmv))]

combat_annot_m = combat_annot_m[,-2]
#match by gene id
combat_annot_m[,c(3:4)] = annot_file[match(combat_annot_m$`colnames(edata_m)`,annot_file$Gene.ID),c(1,2)]

#modNames = substring(names(MEs), 3)
geneModuleMembership_m = as.data.frame(cor(edata_trim_m, MEs_m, use = "p")) # spearman or kendall? using pearson here

MMPvalue_m = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_m), nSamples_m))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership_m$gene = colnames(edata_trim_m)
#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="");


combat_annot_m[5:(ncol(geneModuleMembership_m)+4)] = geneModuleMembership_m[match(geneModuleMembership_m$gene,combat_annot_m$`colnames(edata_m)`),]

x = match(colnames(edata_trim_m),annot_file$Gene.ID)

colnames(edata_trim_m) = annot_file$Gene.Name[x]

#
#
#
#Do the same for females
combat_annot_f = as.data.frame(colnames(edata_f))

combat_annot_f$module = net_f$colors
combat_annot_f$color = moduleColors_f

#REMOVE GREY
rmv = combat_annot_f[which(combat_annot_f$color == "grey"),"colnames(edata_f)"]

combat_annot_f = combat_annot_f[-which(combat_annot_f$`colnames(edata_f)` %in% rmv),]
edata_trim_f = edata_f[,-(which(colnames(edata_f) %in% rmv))]

combat_annot_f = combat_annot_f[,-2]

#change to work for gene id
combat_annot_f[,c(3:4)] = annot_file[match(combat_annot_f$`colnames(edata_f)`,annot_file$Gene.ID),c(1,2)]

#modNames = substring(names(MEs), 3)
geneModuleMembership_f = as.data.frame(cor(edata_trim_f, MEs_f, use = "p")) # spearman or kendall? using pearson here

MMPvalue_f = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership_f), nSamples_f))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership_f$gene = colnames(edata_trim_f)


combat_annot_f[5:(ncol(geneModuleMembership_f)+4)] = geneModuleMembership_f[match(geneModuleMembership_f$gene,combat_annot_f$`colnames(edata_f)`),]

x = match(colnames(edata_trim_f),annot_file$Gene.ID)

colnames(edata_trim_f) = annot_file$Gene.Name[x]

#
#
#
#

save(combat_annot_m, file = "./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
save(combat_annot_f, file = "./results/Rdata/networks/geneModMemAnnot_f_power4.RData")

save(moduleTraitPvalue_f, file = "./results/Rdata/networks/moduleTraitPvalue_f.RData")
save(moduleTraitPvalue_m, file = "./results/Rdata/networks/moduleTraitPvalue_m.RData")

save(moduleTraitCor_f, file = "./results/Rdata/networks/moduleTraitCor_f.RData")
save(moduleTraitCor_m, file = "./results/Rdata/networks/moduleTraitCor_m.RData")

save(edata_trim_m, file = "./results/Rdata/networks/edata_trim_m_5.RData")
save(edata_trim_f, file = "./results/Rdata/networks/edata_trim_f_4.RData")



#########
#########
#Do GO analysis
##
#Females
load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")#combat_annot_f

moduleColors = moduleColors_f[-which(moduleColors_f=="grey")]


network_GO = list()
for(color in unique(combat_annot_f$color)){
  
  allGenes = combat_annot_f$Gene.Name
  interesting.genes = combat_annot_f[which(combat_annot_f$color == color),"Gene.Name"]
  
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
save(network_GO, file="./results/Rdata/networks/GO_Females_sft4.RData")
#####


##
#Males
load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")#combat_annot_m

moduleColors = moduleColors_m[-which(moduleColors_m=="grey")]


network_GO = list()
for(color in unique(combat_annot_m$color)){
  
  allGenes = combat_annot_m$Gene.Name
  interesting.genes = combat_annot_m[which(combat_annot_m$color == color),"Gene.Name"]
  
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
save(network_GO, file="./results/Rdata/networks/GO_Males_sft5.RData")
#####




