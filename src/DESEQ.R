#DESeq on RNA-seq data, for sex and bone strength
library(DESeq2)

#what genes are differentially expressed in hi vs lo bone strength? sex? define hi/lo bone strength by for each sex separately



#create dataset. need to change genes IDs to ENSMUSG or gene names

counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

annot_file = read.csv("~/Documents/projects/DO_project/results/flat/annot_file.csv")
annot_file = annot_file[,c(1,2)]
counts = counts[which(rownames(counts) %in% annot_file$Gene.ID),]

rnames = as.data.frame(rownames(counts))
rnames$gene.name = NA
rnames$gene.name = apply(rnames, 1, function(x) annot_file[which(annot_file$Gene.ID == x[1]),2] ) ##CHECK

rownames(counts) = rnames$gene.name #row=gene, col = ind

##ADD IDENTIFIERS FROM SUPP.R WHEN MAKING RESULTS 
#dont forget to remove low counts. do as before or in deseq manner? probably deseq (actually might not need to)



#######ADD COVARS########
#covars. account for batch, age, sex, and add bone strength hi/lo based on each sex
load("./results/Rdata/cross_eqtl_REDO.Rdata")

#create a covar object from covariates in cross file
#must be numeric
covar = as.matrix(cross_eqtl$covar)
#covar[,"sex"] = (covar[,"sex"] == "M")*1 #convert sex to 1's and 0's
covar[,1] = as.factor(covar[,1]) #sac date to factors
covar[,6] = as.factor(covar[,6]) #generation to factors

#covar = apply(covar,2,as.numeric)
rownames(covar) = rownames(cross_eqtl$covar)
####
rownames(covar)[which(rownames(covar) == "4.1")] = "4"
rownames(covar)[which(rownames(covar) == "371.1")] = "371"

covar = covar[match(colnames(counts),rownames(covar)),]
covar = as.data.frame(covar)
covar$sex = factor(covar$sex)


##ADD BATCH AND CHECK

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


covar$batch1=NA
covar$batch2=NA
covar$batch3=NA
covar$batch4=NA


#define batches
covs$batch1[match(b1,rownames(covs))] = 1
covs[which(rownames(covs)%in% b1==FALSE),"batch1"] = 0

covs$batch2[match(b2,rownames(covs))] = 1
covs[which(rownames(covs)%in% b2==FALSE),"batch2"] = 0

covs$batch3[match(b3,rownames(covs))] = 1
covs[which(rownames(covs)%in% b3==FALSE),"batch3"] = 0

covs$batch4[match(b4,rownames(covs))] = 1
covs[which(rownames(covs)%in% b4==FALSE),"batch4"] = 0

covs$batch = NA

covs[which(covs$batch1 == 1),"batch"] = 1
covs[which(covs$batch2 == 1),"batch"] = 2
covs[which(covs$batch3 == 1),"batch"] = 3
covs[which(covs$batch4 == 1),"batch"] = 4


covs$batch = factor(covs$batch)
covs$sex = factor(covs$sex)
covs$ngen = factor(covs$ngen)
covs$age_factor = factor(covs$age_at_sac_days)

covs$age_bin=NA
covs$age_bin[which(covs$age_at_sac_days <= 85)] = "lo"
covs$age_bin[which(covs$age_at_sac_days > 85)] = "hi"


####add bone strength
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")
pheno = cross_basic$pheno
bone_strength = as.data.frame(pheno[,"bending_max_load"])

rownames(bone_strength)[which(rownames(bone_strength) == "4.1")] = "4"
rownames(bone_strength)[which(rownames(bone_strength) == "371.1")] = "371"
colnames(bone_strength) = "max_load"

#bone_strength = bone_strength[which(rownames(bone_strength) %in% rownames(covs)),]

covs_full = merge(covs, bone_strength, by=0)
rownames(covs_full) = covs_full$Row.names

###split bone strength into hi/lo for each sex
males = covs_full[which(covs_full$sex == 1),]
females = covs_full[which(covs_full$sex == 0),]

median_males = median(males$max_load)
median_females = median(females$max_load)

covs_full$max_load_bin = NA
covs_full$max_load_bin[which( (covs_full$sex == 1)  & (covs_full$max_load <= median_males) )] = "lo"
covs_full$max_load_bin[which( (covs_full$sex == 1)  & (covs_full$max_load > median_males) )] = "hi"

covs_full$max_load_bin[which( (covs_full$sex == 0)  & (covs_full$max_load <= median_females) )] = "lo"
covs_full$max_load_bin[which( (covs_full$sex == 0)  & (covs_full$max_load > median_females) )] = "hi"

#####DO I CREATE SPLIT AMONG SEXES IN BONE STRENGTH BASED ON DISTRIBUTION IN 192 SAMPLES OR OVERALL?
#NORMALIZE BONE STRENGTH? ITS ALREDY PRETTY NORMAL
#IF USING AGE AS NUMERIC, MIGHT NEED TO TRANSFORM, CENTER AND SCALE

#VST/PCA for using this data
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = covs_full, 
                              design = ~sex+batch+age_bin+max_load_bin)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE) #blind so design doesnt matter. seems to be the optimal choice for PCA/viz. Use TRUE if using for downstream applications.
#vsdF <- varianceStabilizingTransformation(dds, blind = FALSE)
#plot in first 2 principal components
# plotPCA(vsdF,intgroup = c("sex"))
# plotPCA(vsd,intgroup = c("batch"))
# plotPCA(vsd,intgroup = c("age_factor"))
# plotPCA(vsd,intgroup = c("age_at_sac_days"))
# plotPCA(vsd,intgroup = c("ngen"))
# plotPCA(vsd,intgroup = c("PEER_factor_1"))


##########factominer#########
rv=rowVars(assay(vsd))
select = order(rv,decreasing = TRUE)[seq_len(min(500,length(rv)))]
dat = t(assay(vsd)[select,])
dat2 = cbind(covs_full,dat)


pca = PCA(dat2[,72:ncol(dat2)],graph=F,scale.unit = F)

fviz_screeplot(pca)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$sex,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="Sex"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$sex,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(3,4),
             legend.title="Sex"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$batch,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="Batch"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$batch,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(3,4),
             legend.title="Batch"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$age_bin,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="Age"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$age_bin,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(3,4),
             legend.title="Age"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$max_load_bin,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="Max Load"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$max_load_bin,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(3,4),
             legend.title="Max Load"
             
)

#DESeq
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = covs_full, 
                              design = ~sex+batch+age_bin+max_load_bin)
DE = DESeq(dds)

res = results(DE)
resOrdered = res[order(res$padj),]

#fold change is hi compared to low
res = results(DE, contrast = c("max_load_bin", "hi","lo"))
resOrdered = res[order(res$padj),]


resLFC_max_load = lfcShrink(DE, type="apeglm",coef = "max_load_bin_lo_vs_hi")

resOrdered_LFC_max_load = resLFC_max_load[order(resLFC_max_load$padj),]

plotMA(resLFC_max_load, ylim=c(-2,2))

x = as.data.frame(resOrdered_LFC_max_load)
x = x[which(x$padj <= 0.05),]

x$symbol = rownames(x)
x$MGI = NA

which(x$symbol %in% MGI_markers$`Marker Symbol` == FALSE) #31 has _isoform in its name, edit manually
x$symbol[31] = "Itgb1bp1"

MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)

x$MGI = apply(x,1,function(z) MGI_markers[which(MGI_markers$`Marker Symbol` == z[6]),1])


########

###
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = covs_full, 
                              design = ~batch+age_bin+sex)
DEsex = DESeq(dds)

res_sex = results(DEsex, contrast = c("sex", "0","1"))
resOrdered = res_sex[order(res_sex$padj),]

resLFC_sex = lfcShrink(DEsex, type="apeglm",coef = "sex_1_vs_0")
resOrdered_LFC_sex = resLFC_sex[order(resLFC_sex$padj),]

# apeglm cite: Zhu, A., Ibrahim, J.G., Love, M.I. (2018) Heavy-tailed prior distributions for
# sequence count data: removing the noise and preserving large differences.
# bioRxiv. https://doi.org/10.1101/303255
plotMA(resLFC_sex, ylim=c(-2,2))


x = as.data.frame(resOrdered_LFC_sex)
x = x[which(x$padj <= 0.05),]

x$symbol = rownames(x)
x$MGI = NA

MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)

#change _isoform to the actual name
x$symbol[grep("isoform", x$symbol)] = sapply(strsplit(x$symbol[grep("isoform", x$symbol)],"_"),"[[",1)

x[which(x$symbol %in% MGI_markers$`Marker Symbol` == FALSE) ,]


x$MGI = apply(x,1,function(z) as.character(MGI_markers[which(MGI_markers$`Marker Symbol` == z[6]),1]))

#checker
z=which(x$MGI == "character(0)")
#

#get MGI ID from synonyms
for(i in which(x$MGI == "character(0)")){
  if(length(MGI_markers[grep(x$symbol[i], x = MGI_markers$`Marker Synonyms (pipe-separated)`),`MGI Accession ID`]) >0){
  x$MGI[i] = MGI_markers[grep(x$symbol[i], x = MGI_markers$`Marker Synonyms (pipe-separated)`),`MGI Accession ID`]
  }
}

x$MGI[z] #CHECK

#add ENSEMBL for AC151602.1
x$ENSEMBL = NA
x[which(x$symbol == "AC151602.1"),"ENSEMBL"] = "ENSMUSG00000118462" #seems to also be LHB
x[which(x$symbol == "AC151602.1"),"MGI"] = "NA" #seems to also be LHB

