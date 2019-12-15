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

#read in VST transformed, quantile-normalized data (produced in ./src/normalize_RNAseq.R)
#load("./results/Rdata/counts_vst_qnorm.Rdata")
#For now try non-quantile normalized data, compare with q-normed

# read in the RNA-seq processed counts file
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#read in an annotation file. This is output from the RNA-seq pipeline
#annot_file = read.delim("./data/314-FarberDO2_S6.gene_abund.tab",header = TRUE) #this looks crappy
#annot_file = readGFF("./results/flat/RNA-seq/mus_stringtie_merged.gtf")
#annot_file = annot_file[which(annot_file$type=="transcript"),]
#remove transcripts not on somatic and X chroms
#chr = c(seq(1:19),"X")
#annot_file = annot_file[which(annot_file$seqid %in% chr),]


#annot_file = annot_file[,c(9,11)]
#annot_file = unique(annot_file)
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


#some genes are duplicated
#counts = t(counts)

#colnames(counts) = annot_file[match(colnames(counts),annot_file$Gene.ID),"Gene.Name"]


#for now, give them a unique ID
#colnames(counts)[which(duplicated(colnames(counts)))] = paste0(colnames(counts)[which(duplicated(colnames(counts)))], "_isoform")

#which(duplicated(colnames(counts)))
#colnames(counts)[3570] = paste0(colnames(counts)[3570],".2")
#colnames(counts)[5342] = paste0(colnames(counts)[5342],".2")
#colnames(counts)[6661] = paste0(colnames(counts)[6661],".3")

#colnames(counts)[4133] = paste0(colnames(counts)[4133],".2")
#colnames(counts)[4134] = paste0(colnames(counts)[4134],".3")
#colnames(counts)[18694] = paste0(colnames(counts)[18694],".2")

#counts = t(counts)


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

covs$male = NA
covs$female = NA
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
modcombat = model.matrix(~as.factor(sex) + as.factor(age_at_sac_days), data=covs)

#batch removal
edata = ComBat(dat=vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)

#transpose the matrix
edata = t(edata)

#quantile normalize

#qnorm by gene
#from https://www.nature.com/articles/nature11401#s1 (FTO genotype BMI Visscher et al)

#for(col in 1:ncol(edata)){
#  edata[,col] = qnorm((rank(edata[,col],na.last="keep")-0.5)/sum(!is.na(edata[,col])))
#}
########################
####
#remove pseudogenes
#edata_sans_ps = edata[,-grep("-ps",colnames(edata))]
#edata_sans_ps = edata_sans_ps[,-grep("-rs",colnames(edata_sans_ps))]
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
sft = pickSoftThreshold(edata, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE,corFnc = "bicor", corOptions = list(maxPOutliers =0.1))
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
#did 4, but try 9
net = blockwiseModules(edata, power = 4,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE,
                       saveTOMs = FALSE,
                       verbose = 3, maxPOutliers = 0.1)
#7 bicor
# net = blockwiseModules(edata, power = 7,
#                        TOMType = "signed", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.15,
#                        numericLabels = TRUE,
#                        saveTOMs = FALSE,
#                        verbose = 3, corType = "bicor", maxPOutliers = 0.1)

## 4: 39 modules not including 0, 6004 genes in module 0
##7: 21 mods no 0, 8814 in mod 0
#saveRDS(net, file="./results/Rdata/networks/wgcna_7_BICOR.RDS")

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

####
####
###
# MEDiss = 1-cor(MEs);
# METree = hclust(as.dist(MEDiss), method = "average");# Plot the result
# sizeGrWindow(7, 6)
# plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
# abline(h=0.15, col = "red")
# 
# # Call an automatic merging function
# merge = mergeCloseModules(exprData = edata, moduleColors, cutHeight = 0.2, verbose = 3)
# # The merged module colors
# newColors = merge$colors;
# # Eigengenes of the new merged modules:
# mergedMEs = merge$newMEs;
# 
# sizeGrWindow(12, 9)
# plotDendroAndColors(net$dendrograms[[1]], colors = cbind(merge$allOK, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
####
####
####

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

save(moduleTraitPvalue, file = "./results/Rdata/networks/moduleTraitPvalue_sexcombined_7_BICOR.RData")
save(moduleTraitCor, file = "./results/Rdata/networks/moduleTraitCor_sexcombined_7_BICOR.RData")

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
MEs = MEs[,-37]
textMatrix = signif(moduleTraitPvalue, 1)
####
####

dim(textMatrix) = dim(moduleTraitCor)
textMatrix = textMatrix[-37,c(34:41)]
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
########################
combat_annot = as.data.frame(colnames(edata))

combat_annot$module = net$colors
combat_annot$color = moduleColors

#REMOVE GREY
rmv = combat_annot[which(combat_annot$color == "grey"),"colnames(edata)"]

combat_annot = combat_annot[-which(combat_annot$`colnames(edata)` %in% rmv),]
edata_trim = edata[,-(which(colnames(edata) %in% rmv))]

#have to run this sometimes?
combat_annot = combat_annot[,-2]
#the gsub allows for matching of genes that had _isoform* added to them
#combat_annot[,c(3:4)] = annot_file[match(gsub(combat_annot$`colnames(edata)`,pattern = "_isoform.*",replacement = ""),annot_file$gene_name),c(1,2)]
combat_annot[,c(3:4)] = annot_file[match(combat_annot$`colnames(edata)`,annot_file$Gene.ID),c(1,2)]


#modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(edata_trim, MEs, use = "p")) # spearman or kendall? using pearson here

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership$gene = colnames(edata_trim)
#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="");


combat_annot[5:(ncol(geneModuleMembership)+4)] = geneModuleMembership[match(geneModuleMembership$gene,combat_annot$`colnames(edata)`),]


save(combat_annot, file = "./results/Rdata/networks/geneModMemAnnot_sexcombined_power7_BICOR.RData")

save(moduleTraitPvalue, file = "./results/Rdata/networks/moduleTraitPvalue_sexcombined_power7_BICOR.RData")

save(moduleTraitCor, file = "./results/Rdata/networks/moduleTraitCor_sex_combined_power7_BICOR.RData")







##
load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
combat_annot = geneModMemAnnot
#annot_file_resid[which(moduleColors=="red"),]
moduleColors = moduleColors[-which(moduleColors=="grey")] #check that this didnt screw up the hub calculation on following line
hubs = chooseTopHubInEachModule(edata,moduleColors)
####topGO###
#allGenes = select(Mus.musculus, "ENSMUSG00000066438", columns, "ENSEMBL")
#allGenes = allGenes$ENSEMBL
#interesting.genes<-res$ENSEMBL
#allGenes = colnames(edata)


network_GO = list()
for(color in unique(combat_annot$color)){
  
  allGenes = combat_annot$gene
  interesting.genes = combat_annot[which(combat_annot$color == color),"gene"]
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
  ######
  network_GO[[color]] = t.all
}
save(network_GO, file="./results/Rdata/networks/GO_male_5.RData")
#####
#####
####
allGenes = combat_annot$gene
interesting.genes = combat_annot[which(combat_annot$color == "salmon"),"gene"]
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



sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.35,
               zlim = c(-1,1),
               cex.lab.x = 0.7,
               cex.lab.y = 0.8,
               verticalSeparator.x = c(1:length(names(datTraits))),
               horizontalSeparator.y = c(1:length(names(MEs))),
               main = paste("Module-trait relationships"))
###
###
geneModMemAnnot = combat_annot

#write out net, edata connect and annot. 
#For blacklist, we would use cis-eqtl
save(edata_trim, file = "./results/Rdata/edata_7_BICOR.RData")
save(geneModMemAnnot, file = "./results/Rdata/geneModMemAnnot_power7_BICOR.RData")

########################## CONSTRUCT BAYESIAN NETWORKS FOR EACH MODULE ###########################
#Done on Rivanna. learn_bn.R. Only need to know number of colors (modules) for SLURM script

resid = edata
# #module membership for all genes in all modules (from auto_wgcna.R)
 connect = geneModMemAnnot
# 
 clr = unique(connect$color)
# 
# mod_genes = connect[which(connect$color==i),"gene"]
# 
 #mod_genes_membership = connect[which(connect$`colnames(edata)` %in% mod_genes),c("gene","Gene.ID",paste0("ME",i),"color")]
 #mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == i),]
# 
 #mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$gene)])
# 
 #bn = mmhc(mod_genes_exp)
# 
 #varName = paste0("hybrid_",i,"_nobl")
# 
 #assign(x = varName, value = bn)
 #print(i)
 for(i in clr){
   mod_genes = connect[which(connect$color==i),"gene"]
#   
   mod_genes_membership = connect[which(connect$`colnames(edata)` %in% mod_genes),c("gene","Gene.ID",paste0("ME",i),"color")]
   mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == i),]
#   
   mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$gene)])
# 
   bn = mmhc(mod_genes_exp)
#   
   varName = paste0("hybrid_",i,"_nobl")
#   
   assign(x = varName, value = bn)
   print(i)
 }
# 
# 
# 
# 
# 
# 
# 
# 
# #PHENOTYPE DATA
# #fix to read object in directly or link to file
# # full_pheno_table_ordered = read.delim("~/Desktop/DO_proj/pheno_data/full_pheno_table_ordered",stringsAsFactors = FALSE,as.is = TRUE)
# # datTraits = full_pheno_table_ordered[which(rownames(full_pheno_table_ordered) %in% rownames(resid)),]
# # #datTraits[,135] = as.numeric(gsub(x = datTraits[,135],pattern = "~ ",replacement = ""))
# # datTraits = datTraits[,-c(1:7,16,72:150)]
# # for(i in 1:ncol(datTraits)){datTraits[,i] = as.numeric(datTraits[,i])}
# #datTraits = log2(datTraits + 0.0001)
# 
# 
# 
# 
# ############## Assign blacklist based on eQTL #############
# ############## Mouse GWS is 5*10^-5 #######################
# #the commented out portion used old microarray data. Use DO cis eqtl now from (cis_and_trans_eqtl.R)
# # cis_eSNPs = read.csv("cis.eSNP.031609.csv", as.is = TRUE, stringsAsFactors = FALSE)
# # cis_eSNPs_ME9 = cis_eSNPs[which(cis_eSNPs$Probe %in% m9_genes_most_connected$TargetId),]
# # cis_eSNPs_ME9$bp = as.numeric(cis_eSNPs_ME9$bp)*1000000
# # 
# # #some gene names are different between cis_eSNPs_ME9 and m9_unique colnames
# # cis_eSNPs_ME9_gws = cis_eSNPs_ME9[which(cis_eSNPs_ME9$cisp <= 5e-5),]
# 
# cis_genes = unique(cis_eqtl_ALL$lodcolumn)
# 
# #blacklist_genes_parents = m9_genes_most_connected[which(m9_genes_most_connected$TargetId %in% cis_eSNPs_ME9_gws$Probe),"Symbol"]
# blacklist_genes_parents = mod_genes[which(mod_genes %in% cis_genes)]
# 
# #a gene not in blacklist_genes cannot be a parent to a gene in blacklist_genes
# #blacklist_genes_children = colnames(m9_genes_exp)[which(colnames(m9_genes_exp) %in% blacklist_genes_parents == F)]
# blacklist_genes_children = mod_genes[which(mod_genes %in% cis_genes == FALSE)]
# 
# blacklist_df = as.data.frame(rep(blacklist_genes_children,times = length(blacklist_genes_parents)),stringsAsFactors = F)
# colnames(blacklist_df)[1] = "From"
# blacklist_df$To = rep(blacklist_genes_parents, each = length(blacklist_genes_children))
# 
# #pheno_bl = as.data.frame(rep("mean.bmd.sp",times = 354),stringsAsFactors = F)
# #colnames(pheno_bl)[1] = "From"
# #pheno_bl$To = rep(colnames(m9_exp_pheno[,c(1:354)]), each = 1)
# 
# #blacklist_df = rbind(blacklist_df, pheno_bl)
# ############################################################
# ############################################################
# ############################################################
# #red_genes_exp[] = lapply(red_genes_exp, as.double)
# #hybrid_new_test = mmhc(m9_exp_pheno[,c(1:354,356)],blacklist = blacklist_df)
# #hybrid_new_test_nobl = mmhc(m9_exp_pheno[,c(1:354,356)])
# 
# 
# #remove duplicated genes (same module membership, for now just removing ones that appear second. Different expression patterns though)
# #x = which(colnames(brown_genes_exp) %in% colnames(brown_genes_exp)[which(duplicated(colnames(brown_genes_exp)))])
# #x = brown_genes_exp[,x]
# 
# #test if duplicates exist. there shouldnt be any
# which(duplicated(colnames(mod_genes_exp)))
# #remove if dups exist (or rename)
# #mod_genes_exp = mod_genes_exp[,-which(duplicated(colnames(mod_genes_exp)))]
# 
# #colnames(mod_genes_exp)[2357] = "Wscd2.dup"
# colnames(mod_genes_exp)[993] = "Mia3.dup2"
# hybrid_blue_nobl = mmhc(mod_genes_exp)#better BIC (higher values better, AIC and BIC rescaled by 2)
# bnlearn::score(hybrid_blue_nobl, mod_genes_exp)
# 
# hybrid_blue_bl = mmhc(mod_genes_exp,blacklist = blacklist_df)
# score(hybrid_blue_bl, mod_genes_exp)
# 
# 
mod_genes_membership = geneModMemAnnot[which(geneModMemAnnot$color=="brown"),]
mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$gene)])
#colnames(mod_genes_exp)[which(colnames(mod_genes_exp)=="Ppp2r3d_isoform")] = "Ppp2r3d_isoform"
hybrid_arcs = arc.strength(bn, data = mod_genes_exp)
hybrid_arcs = strength.plot(hybrid_new_test, hybrid_arcs)
boot_arcs = boot.strength(data = mod_genes_exp,cpdag=T, algorithm = "mmhc" )
attr(boot_arcs, "threshold")
ave = averaged.network(boot_arcs)


#plot MEred with Tb.N

tbn = datTraits$uCT_Tb.N
red = MEs$MEred
x = as.data.frame(t(rbind(tbn,red)))

ggplot(x, aes(x = tbn, y = red)) +
  geom_point(aes(col = "#C42126")) +theme_bw() +
  theme(plot.title = element_text(size = 12, family = "Tahoma", face = "bold", hjust=0.5),
        text = element_text(size = 12, family = "Tahoma"), legend.position = "none") +
  scale_fill_brewer(palette="Dark2") + labs(title = "Red Module Eigengene vs. Trabecular Number",x="Trabecular Number", y="MERed Eigengene")+geom_smooth(method = "lm", color="red")



