############
library(sva)
library(Mus.musculus)#seems to use mm10
library(WGCNA)
library(topGO)
library(DESeq2)
library("igraph")
library("bnlearn")
library("parallel")
############
options(stringsAsFactors = FALSE)
#

#read in VST transformed, quantile-normalized data (produced in ./src/normalize_RNAseq.R)
#load("./results/Rdata/counts_vst_qnorm.Rdata")
#For now try non-quantile normalized data, compare with q-normed

# read in the RNA-seq processed counts file
counts = read.csv("./results/flat/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#read in an annotation file. This is output from the RNA-seq pipeline
annot_file = read.delim("./data/314-FarberDO2_S6.gene_abund.tab",header = TRUE)
#
#remove transcripts not on somatic and X chroms
chr = c(seq(1:19),"X")
rmv = annot_file[which(annot_file$Reference %in% chr == FALSE),"Gene.ID"]

counts = counts[-which(rownames(counts) %in% rmv),]
#

#find and remove features that have fewer than 10 reads in more than 90% (173) of samples 
x=c()
for(i in 1:nrow(counts)){
  if(sum(counts[i,]<10)>=173){
    print(i)
    x = append(x,i)
  }
}

#253 genes removed
counts = counts[-x,]


#some genes are duplicated
counts = t(counts)

colnames(counts) = annot_file[match(colnames(counts),annot_file$Gene.ID),"Gene.Name"]
counts[,colnames(counts)[which(duplicated(colnames(counts)))]]

#for now, give them a unique ID
colnames(counts)[which(duplicated(colnames(counts)))] = paste0(colnames(counts)[which(duplicated(colnames(counts)))], "_isoform")

which(duplicated(colnames(counts)))
colnames(counts)[4131] = paste0(colnames(counts)[4131],".2")
colnames(counts)[4132] = paste0(colnames(counts)[4132],".3")
colnames(counts)[18687] = paste0(colnames(counts)[18687],".2")

counts = t(counts)


######


#vst from deseq2
vst = varianceStabilizingTransformation(as.matrix(counts))

#get batch. What I did here is I got batch from the file names of the alignment output for RNA-seq
f = list.files("./results/flat/RNA-seq/sums/")
p1 = f[grep(pattern = "Pool1",f)]
p2 = f[grep(pattern = "Pool2",f)]
p3 = f[grep(pattern = "FarberDO2",f)]

b1 = c()
b2 = c()
b3 = c()
for(i in 1:length(p1)){
  b1 = c(b1,strsplit(strsplit(p1[i],"Pool1_")[[1]][2],"-")[[1]][1])
  b1 = unique(b1)
}

for(i in 1:length(p2)){
  b2 = c(b2,strsplit(strsplit(p2[i],"Pool2_")[[1]][2],"-")[[1]][1])
  b2 = unique(b2)
}

for(i in 1:length(p3)){
  b3 = c(b3,strsplit(strsplit(p3[i],"fastq_")[[1]][2],"-")[[1]][1])
  b3 = unique(b3)
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


#covs$male = as.numeric(covs$sex=="M")
#covs$female = as.numeric(covs$sex=="F")
#define batches
covs$batch1[match(b1,rownames(covs))] = 1
covs[which(rownames(covs)%in% b1==FALSE),"batch1"] = 0

covs$batch2[match(b2,rownames(covs))] = 1
covs[which(rownames(covs)%in% b2==FALSE),"batch2"] = 0

covs$batch3[match(b3,rownames(covs))] = 1
covs[which(rownames(covs)%in% b3==FALSE),"batch3"] = 0


#cc = as.matrix(covs[,c(3,5:8)])
#rownames(cc) = rownames(covs)

#convert to integers
#covs[,1] = as.factor(covs[,1])

#cc[,c(1:5)] = as.integer(cc[,])

covs = as.data.frame(covs[,2:9])

covs$batch = NA
covs[which(covs$batch1 == 1),"batch"] = 1
covs[which(covs$batch2 == 1),"batch"] = 2
covs[which(covs$batch3 == 1),"batch"] = 3


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
#did 4, but try 9
net = blockwiseModules(edata, power = 9,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)

## 35 modules not including 0, 5736 genes in module 0
#for thresh of 9, 20 mod, 10572 in mod 0
#thresh 6, 6717 in mod 0, 31 modules
#
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
#cor module eigengenes with traits
moduleTraitCor = cor(MEs, datTraits, use = "p",method = "p");

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


moduleTraitPvalue = as.matrix(moduleTraitPvalue)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

trait_names = gsub(pattern = "bending_",replacement = "",x = names(datTraits))
trait_names = gsub(pattern = "uCT_",replacement = "",x = trait_names)
trait_names = gsub(pattern = "\\.\\.",replacement = "\\.",x = trait_names)
trait_names = gsub(pattern = "\\.\\.",replacement = "\\.",x = trait_names)
trait_names[8] = "Adiposity"


# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = trait_names,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               zlim = c(-1,1),
               cex.lab.x = 0.7,
               cex.lab.y = 0.8,
               verticalSeparator.x = c(1:length(names(datTraits))),
               horizontalSeparator.y = c(1:length(names(MEs))),
               main = paste("Module-trait relationships"))
###
which(moduleTraitCor == sort((moduleTraitCor),decreasing =TRUE)[1], arr.ind = T)
colnames(moduleTraitCor)[33]
########################
combat_annot = as.data.frame(colnames(edata))

combat_annot$module = net$colors
combat_annot$color = moduleColors

#the gsub allows for matching of genes that had _isoform* added to them
combat_annot[,c(4:8)] = annot_file[match(gsub(combat_annot$`colnames(edata)`,pattern = "_isoform.*",replacement = ""),annot_file$Gene.Name),c(1,3,4,5,6)]


modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(edata, MEs, use = "p"))

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership$gene = colnames(edata)
#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="");

combat_annot[9:(ncol(geneModuleMembership)+8)] = geneModuleMembership[match(geneModuleMembership$gene,combat_annot$`colnames(edata)`),c(1:ncol(geneModuleMembership))]


#which(moduleColors=="red")
#annot_file_resid[which(moduleColors=="red"),]
hubs = chooseTopHubInEachModule(edata,moduleColors)
####topGO###
#allGenes = select(Mus.musculus, "ENSMUSG00000066438", columns, "ENSEMBL")
#allGenes = allGenes$ENSEMBL
#interesting.genes<-res$ENSEMBL
#allGenes = colnames(edata)
allGenes = combat_annot$gene
interesting.genes = combat_annot[which(combat_annot$color == "brown"),"gene"]
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
save(edata, file = "./results/Rdata/edata.RData")
save(geneModMemAnnot, file = "./results/Rdata/geneModMemAnnot_power9.RData")

########################## CONSTRUCT BAYESIAN NETWORKS FOR EACH MODULE ###########################
#Done on Rivanna. learn_bn.R. Only need to know number of colors (modules) for SLURM script

# resid = edata
# #module membership for all genes in all modules (from auto_wgcna.R)
# connect = geneModMemAnnot
# 
# clr = unique(connect$color)
# 
# mod_genes = connect[which(connect$color==i),"gene"]
# 
# mod_genes_membership = connect[which(connect$`colnames(edata)` %in% mod_genes),c("gene","Gene.ID",paste0("ME",i),"color")]
# mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == i),]
# 
# mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$gene)])
# 
# bn = mmhc(mod_genes_exp)
# 
# varName = paste0("hybrid_",i,"_nobl")
# 
# assign(x = varName, value = bn)
# print(i)
# for(i in clr){
#   mod_genes = connect[which(connect$color==i),"gene"]
#   
#   mod_genes_membership = connect[which(connect$`colnames(edata)` %in% mod_genes),c("gene","Gene.ID",paste0("ME",i),"color")]
#   mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == i),]
#   
#   mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$gene)])
# 
#   bn = mmhc(mod_genes_exp)
#   
#   varName = paste0("hybrid_",i,"_nobl")
#   
#   assign(x = varName, value = bn)
#   print(i)
# }
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
# hybrid_arcs = arc.strength(hybrid_blue_nobl, data = mod_genes_exp)
# #strength.plot(hybrid_new_test, hybrid_arcs)
