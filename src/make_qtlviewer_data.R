library(qtl2convert)
library(qtl2)
library(tidyverse)
options(scipen=999,stringsAsFactors = FALSE,digits = 15)
#make Rdata file for qtlviewer

# The following elements should be contained within the RData file.
# ensembl.version – the numerical version of Ensembl
# genoprobs – the genotype probabilities
# K – the kinship matrix
# map – list of one element per chromosome, with the genomic position of each marker
# markers – marker names and positions
# The following element is a special element and there must be at least one per RData file.
# dataset.* - where * should be a very short, unique and informative name. This element
# will contain most of the data and will be detailed in the section below.
# Exact case of element and variable names is very important.
# Other meta data can be included in the RData file as long as there are no conflicting names.

ensembl.version = 94 #check this (97)
#
load(file = "./results/Rdata/apr_basic_cleaned.Rdata") #genoprobs
genoprobs = apr
rm(apr)
#
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
K = k_loco
rm(k_loco)
#
load(file = "./results/Rdata/cross_eqtl.Rdata")
map = cross_eqtl$pmap
#
markers = qtl2convert::map_list_to_df(cross_eqtl$pmap)
colnames(markers)[3] = "marker.id"
markers = as_tibble(markers)

##########
##########
# The environment must contain at least one object of this type, multiple are allowed. The *
#   should be a very short, unique and informative name. It is for internal use only and will not
# appear in the QTL Viewer interface.
# The main purpose of the dataset.* element is to store multiple datasets per RData file with
# informative information regarding the data.

# The dataset.* element is a list that should contain the following named elements,
# annot.datatype – annotations, where datatype is one of mrna, protein, or phenotype
# annot.samples – annotation data for the samples
# covar.matrix – a matrix of covariate data, samples (rows) x covariates (columns)
# covar.info – information describing the covariates
# data – either a matrix containing data or a list containing several kinds of data
# datatype – one of mrna, protein, or phenotype
# display.name – name of the dataset, for QTL Viewer display purposes
# lod.peaks – a list of LOD peaks over a certain threshold

#cis
gene_annot_file = read.table("~/Desktop/DO_proj/data/GIGAMUGA/gene_annot_file", header = TRUE,stringsAsFactors = FALSE)
annot.mrna = gene_annot_file
#colnames(annot.mrna) = c("gene.id","symbol","chr","start","TSS","end")
colnames(annot.mrna) = c("gene.id","symbol","chr","start","end")

annot.mrna$start = annot.mrna$start/1000000
annot.mrna$end = annot.mrna$end/1000000
#annot.mrna$TSS = annot.mrna$TSS/1000000

chr = c(1:19, "X")
annot.mrna[which(annot.mrna$chr %in% chr == FALSE),] = NA
annot.mrna = na.omit(annot.mrna)
annot.mrna$middle = (floor((annot.mrna$start*1000000+annot.mrna$end*1000000)/2))/1000000
annot.mrna$nearest.marker.id = NA
for(i in 1:nrow(annot.mrna)){
  
  annot.mrna$nearest.marker.id[i] = find_marker(map = map,chr = annot.mrna$chr[i],pos = annot.mrna$middle[i])
}
#
annot.samples = cross_eqtl$covar
annot.samples$mouse.id = rownames(annot.samples)
annot.samples = annot.samples[c(65,1:64)]
annot.samples$sex = (annot.samples$sex == "M")*1
#
covar.matrix = cross_eqtl$covar
covar.matrix = covar.matrix[,c(2,17:51)]#only sex and first 35
covar.matrix$sex = (covar.matrix$sex == "M")*1
colnames(covar.matrix)[1] = "sexM"
covar.matrix = as.data.frame(lapply(covar.matrix, as.numeric))
covar.matrix = as.matrix(covar.matrix)
rownames(covar.matrix) = rownames(cross_eqtl$covar)
#
covar.info = as.data.frame(matrix(ncol=6,nrow = 2))
colnames(covar.info) = c("sample.column","covar.column","display.name","interactive","primary","lod.peaks")
covar.info[1,] = c("sex","sexM","Sex",FALSE,FALSE,NA)
covar.info[2,] = c("ngen",NA,"Generation",FALSE,FALSE,NA)
#
annot = read.delim("~/Desktop/DO_proj/data/networks/314-FarberDO2_S6.gene_abund.tab",stringsAsFactors = FALSE)

data = cross_eqtl$pheno #must be ensembl ids
gene_names = annot[match(colnames(data),annot$Gene.ID),"Gene.Name"]
#gene_names = gene_annot_file[match(colnames(data),gene_annot_file$mgi_symbol),"ensembl_gene_id"]

#this adds 
ensembl_names = gene_annot_file[match(gene_names,gene_annot_file$mgi_symbol),"ensembl_gene_id"]

colnames(data) = ensembl_names#not all unique

data = data[,-which(colnames(data) %in% annot.mrna$gene.id == FALSE)]

data = data[,-which(duplicated(colnames(data)))]
#this adds colnames that arent in ensembl annotaiton. dont use because need 1:1 mapping
#colnames(data)[which(is.na(colnames(data)))] = gene_names[which(is.na(colnames(data)))]

data = as.matrix(data)
#
datatype = "mrna"
display.name = "Local eQTL"
#
lod.peaks = c()
lod.peaks = as.list(lod.peaks)
lod.peaks$additive = NA

load("./results/Rdata/cis_eqtl_ALL.Rdata")

lod.peaks$additive = cis_eqtl_ALL

names(lod.peaks$additive)[8] = "gene.id"
lod.peaks$additive$marker.id = NA

for(i in 1:nrow(lod.peaks$additive)){
  lod.peaks$additive$marker.id[i] = find_marker(map = map,chr = lod.peaks$additive$chr[i],pos = lod.peaks$additive$pos[i])
  
}

dataset.localEqtl = list()

dataset.localEqtl[[1]] = annot.mrna
dataset.localEqtl[[2]] = annot.samples
dataset.localEqtl[[3]] = covar.matrix
dataset.localEqtl[[4]] = covar.info
dataset.localEqtl[[5]] = data
dataset.localEqtl[[6]] = datatype
dataset.localEqtl[[7]] = display.name
dataset.localEqtl[[8]] = lod.peaks

names(dataset.localEqtl) = c("annot.mrna", "annot.samples","covar.matrix","covar.info","data","datatype","display.name","lod.peaks")


###do same for distalEqtl. Difference will be that covar is different (all and no sex) and lod peaks will be different. Also display name
###
###
###

display.name = "Distal eQTL"

covar.matrix = cross_eqtl$covar
covar.matrix = covar.matrix[,c(17:64)]#all PEER, no sex

covar.matrix = as.data.frame(lapply(covar.matrix, as.numeric))
covar.matrix = as.matrix(covar.matrix)
rownames(covar.matrix) = rownames(cross_eqtl$covar)
#
covar.info = as.data.frame(matrix(ncol=6,nrow = 2))
colnames(covar.info) = c("sample.column","covar.column","display.name","interactive","primary","lod.peaks")
covar.info[1,] = c("sex",NA,"Sex",FALSE,FALSE,NA)
covar.info[2,] = c("ngen",NA,"Generation",FALSE,FALSE,NA)
#


lod.peaks = c()
lod.peaks = as.list(lod.peaks)
lod.peaks$additive = NA

load("./results/Rdata/trans_eqtl_ALL.Rdata")

lod.peaks$additive = trans_eqtl_ALL

names(lod.peaks$additive)[8] = "gene.id"
lod.peaks$additive$marker.id = NA

for(i in 1:nrow(lod.peaks$additive)){
  lod.peaks$additive$marker.id[i] = find_marker(map = map,chr = lod.peaks$additive$chr[i],pos = lod.peaks$additive$pos[i])
  
}



#
dataset.distalEqtl = list()

dataset.distalEqtl[[1]] = annot.mrna
dataset.distalEqtl[[2]] = annot.samples
dataset.distalEqtl[[3]] = covar.matrix
dataset.distalEqtl[[4]] = covar.info
dataset.distalEqtl[[5]] = data
dataset.distalEqtl[[6]] = datatype
dataset.distalEqtl[[7]] = display.name
dataset.distalEqtl[[8]] = lod.peaks

names(dataset.distalEqtl) = c("annot.mrna", "annot.samples","covar.matrix","covar.info","data","datatype","display.name","lod.peaks")

#

# rm(annot_mstrg)
# rm(gene_annot_file)
# rm(peaks_35peer_qnorm_ALL)
# rm(x1,x2,x3,x4)
# rm(y1,y2,y3,y4)
# rm(Xcovar)
# rm(getCisEqtl)
# rm(getTransEqtl)
# rm(gene_names)
# rm(i)
# rm(chr)
# rm(cis_eqtl_ALL)
# rm(trans_eqtl_ALL)
# rm(ensembl_names)
# rm(annot,annot.mrna,annot.samples,covar.info,covar.matrix,cross_eqtl,data)
# rm(lod.peaks)
# rm(covar)
# rm(peaks_ALLpeer_qnorm_ALL)
#save.image("~/Desktop/example/test.RData")

#pheno
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")
map = cross_basic$pmap
#
#load("~/Desktop/example.RData")
#
annot.phenotype = as.data.frame(matrix(nrow = ncol(cross_basic$pheno)+10,ncol = 17))
colnames(annot.phenotype) = c("data.name","short.name","R.name","description","units","is.id","category","R.category","is.numeric","is.date","is.factor","factor.levels","is.covar","is.pheno","is.derived","omit","use.covar")
annot.phenotype$data.name = c(colnames(cross_basic$pheno),"generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")
annot.phenotype$data.name[1] = "mouse.id"

annot.phenotype$short.name = c(colnames(cross_basic$pheno),"generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")
annot.phenotype$short.name[1] = "mouse.id"

annot.phenotype$R.name = c(colnames(cross_basic$pheno),"generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")
annot.phenotype$R.name[1] = "mouse.id"

annot.phenotype$description= c(colnames(cross_basic$pheno),"generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")
annot.phenotype$units = NA
annot.phenotype$is.id= c(TRUE, rep(FALSE,86))
annot.phenotype$category = NA
annot.phenotype$R.category = NA
annot.phenotype$is.numeric = rep(TRUE,87)
annot.phenotype$is.date = rep(FALSE,87)
annot.phenotype$is.factor = rep(FALSE,87)
annot.phenotype$factor.levels = rep(NA,87)
annot.phenotype$is.covar = c(NA,NA,TRUE,TRUE,TRUE,rep(NA,72),rep(TRUE,10))
annot.phenotype$is.pheno = c(FALSE,FALSE,FALSE,FALSE,FALSE,rep(TRUE,72),rep(FALSE,10))
annot.phenotype$is.derived = rep(FALSE,87)
annot.phenotype$omit = rep(FALSE,87)
annot.phenotype$use.covar = c(rep(NA,5),rep("sex:age_at_sac_days:body_weight:generationG24:generationG25:generationG26:generationG27:generationG28:generationG29:generationG30:generationG31:generationG32:generationG33",72),rep(NA,10))
####
annot.samples = as.data.frame(matrix(nrow = nrow(cross_basic$pheno),ncol = 15))
colnames(annot.samples) = c("mouse.id","sex","sac_date","age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")

covar = as.data.frame(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1
covar[,2:16] = apply(covar[,2:16],2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file

annot.samples$mouse.id = rownames(covar)
annot.samples$sex = covar$sex
annot.samples$sac_date = covar$sac_date
annot.samples$age_at_sac_days = covar$age_at_sac_days
annot.samples$body_weight = covar$body_weight
annot.samples$generationG24 = covar$generationG24
annot.samples$generationG25 = covar$generationG25
annot.samples$generationG26 = covar$generationG26
annot.samples$generationG27 = covar$generationG27
annot.samples$generationG28 = covar$generationG28
annot.samples$generationG29 = covar$generationG29
annot.samples$generationG30 = covar$generationG30
annot.samples$generationG31 = covar$generationG31
annot.samples$generationG32 = covar$generationG32
annot.samples$generationG33 = covar$generationG33

##########
display.name = "Phenotypic QTL"

covar.matrix = covar[c(2,3,4,7:16)]

covar.matrix = as.data.frame(lapply(covar.matrix, as.numeric))
covar.matrix = as.matrix(covar.matrix)
rownames(covar.matrix) = rownames(covar)
covar.matrix = as.matrix(na.omit(covar.matrix))

#
covar.info = as.data.frame(matrix(ncol=6,nrow = ncol(covar.matrix)))
colnames(covar.info) = c("sample.column","covar.column","display.name","interactive","primary","lod.peaks")
covar.info[1,] = c("sex","sex","Sex",TRUE,TRUE,"sexInt")
covar.info[2,] = c("age_at_sac_days","age_at_sac_days","Age_in_days",FALSE,TRUE,NA)
covar.info[3,] = c("body_weight","body_weight","Body_weight",FALSE,TRUE,NA)
covar.info[4,] = c("generationG24","generationG24","Generation24",FALSE,TRUE,NA)
covar.info[5,] = c("generationG25","generationG25","Generation25",FALSE,TRUE,NA)
covar.info[6,] = c("generationG26","generationG26","Generation26",FALSE,TRUE,NA)
covar.info[7,] = c("generationG27","generationG27","Generation27",FALSE,TRUE,NA)
covar.info[8,] = c("generationG28","generationG28","Generation28",FALSE,TRUE,NA)
covar.info[9,] = c("generationG29","generationG29","Generation29",FALSE,TRUE,NA)
covar.info[10,] = c("generationG30","generationG30","Generation30",FALSE,TRUE,NA)
covar.info[11,] = c("generationG31","generationG31","Generation31",FALSE,TRUE,NA)
covar.info[12,] = c("generationG32","generationG32","Generation32",FALSE,TRUE,NA)
covar.info[13,] = c("generationG33","generationG33","Generation33",FALSE,TRUE,NA)

####
raw = cross_basic$pheno
raw = cbind(raw,covar[,7:16])
raw$sex = covar$sex

norm_pheno = as.data.frame(log10(cross_basic$pheno[,c(6:14,16,17,21,23:32,34,35,37:41,43:49,51,52,54:58,60:70,72,74,76)]))
pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,33,36,42,50,53,59)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}

log = pheno_combined
log = cbind(log,raw[,colnames(raw)[which(colnames(raw) %in% colnames(pheno_combined)==FALSE)]])
log = log[,colnames(raw)]
log$sex = covar$sex
##remove non log normalized values in log frame to remove confusion? would laso have to remove them from the LOD score frame

#data = list()
#data[[1]] = raw
#data[[2]] = log
data = log
data = as.matrix(data)
#names(data) = c("raw","log")
##
datatype = "phenotype"
##
#lod.peaks is a list with "additive" and "sexInt" elements
load("./results/Rdata/DO_qtl_scan_norm.Rdata")
qtl_peaks_norm = find_peaks(DO_qtl_scan_normal, cross_basic$pmap, threshold=4, drop=1.5)

load("./results/Rdata/DO_qtl_scan_bin_norm.Rdata")
qtl_peaks_bin_norm = find_peaks(DO_qtl_scan_binary_norm, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_both_norm = rbind(qtl_peaks_norm,qtl_peaks_bin_norm)


additive = qtl_peaks_both_norm
  
  for(i in 1:nrow(additive)){
    additive$marker.id[i] = find_marker(map = map,chr = additive$chr[i],pos = additive$pos[i])
    
  }
additive = additive[,c(2,8,5)]
colnames(additive) = c("data.name","marker.id","lod")
#

load("./results/Rdata/DO_qtl_int.Rdata")
qtl_peaks_int = find_peaks(DO_qtl_int, cross_basic$pmap, threshold=4, drop=1.5)

load("./results/Rdata/DO_qtl_int_binary.Rdata")
qtl_peaks_binary_int = find_peaks(DO_qtl_int_binary, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_both_int = rbind(qtl_peaks_int,qtl_peaks_binary_int)

sexInt = qtl_peaks_both_int

for(i in 1:nrow(sexInt)){
  sexInt$marker.id[i] = find_marker(map = map,chr = sexInt$chr[i],pos = sexInt$pos[i])
  
}
sexInt = sexInt[,c(2,8,5)]
colnames(sexInt) = c("data.name","marker.id","lod")
#
lod.peaks = list()
lod.peaks[[1]] = additive
lod.peaks[[2]] = sexInt
names(lod.peaks) = c("additive","sexInt")

########


dataset.phenotype = list()

dataset.phenotype[[1]] = annot.phenotype
dataset.phenotype[[2]] = annot.samples
dataset.phenotype[[3]] = covar.matrix
dataset.phenotype[[4]] = covar.info
dataset.phenotype[[5]] = data
dataset.phenotype[[6]] = datatype
dataset.phenotype[[7]] = display.name
dataset.phenotype[[8]] = lod.peaks

names(dataset.phenotype) = c("annot.phenotype", "annot.samples","covar.matrix","covar.info","data","datatype","display.name","lod.peaks")







###convert to tiblle (from MATT)

temp <- list(
  
  annot.mrna      = tibble::as_tibble(dataset.distalEqtl$annot.mrna),
  
  annot.samples   = tibble::as_tibble(dataset.distalEqtl$annot.samples),
  
  covar.matrix    = dataset.distalEqtl$covar.matrix,
  
  covar.info      = tibble::tibble(
    
    sample.column = c("sex", "ngen"),
    
    display.name = c("Sex", "Generation"),
    
    interactive = c(FALSE, FALSE),
    
    primary = c(TRUE, FALSE),
    
    lod.peaks = c(NA, NA),
    
  ),
  
  data            = dataset.distalEqtl$data,
  
  datatype        = "mRNA",
  
  display.name    = dataset.distalEqtl$display.name,
  
  lod.peaks       = list(additive = tibble::as_tibble(dataset.distalEqtl$lod.peaks$additive))
  
)

dataset.distalEqtl <- temp



temp <- list(
  
  annot.mrna      = tibble::as_tibble(dataset.localEqtl$annot.mrna),
  
  annot.samples   = tibble::as_tibble(dataset.localEqtl$annot.samples),
  
  covar.matrix    = dataset.localEqtl$covar.matrix,
  
  covar.info      = tibble::tibble(
    
    sample.column = c("sex", "ngen"),
    
    display.name = c("Sex", "Generation"),
    
    interactive = c(FALSE, FALSE),
    
    primary = c(TRUE, FALSE),
    
    lod.peaks = c(NA, NA),
    
  ),
  
  data            = dataset.localEqtl$data,
  
  datatype        = "mRNA",
  
  display.name    = dataset.localEqtl$display.name,
  
  lod.peaks       = list(additive = tibble::as_tibble(dataset.localEqtl$lod.peaks$additive))
  
)

dataset.localEqtl <- temp







temp <- list(
  
  annot.phenotype = tibble::as_tibble(dataset.phenotype$annot.phenotype),
  
  annot.samples   = tibble::as_tibble(dataset.phenotype$annot.samples),
  
  covar.matrix    = dataset.phenotype$covar.matrix,
  
  covar.info      = dataset.phenotype$covar.info,
  
  #data            = list(raw = as.matrix(dataset.phenotype$data$raw),
                         
   #                      log = as.matrix(dataset.phenotype$data$log)),
  
  data = dataset.phenotype$data,
  
  datatype        = "phenotype",
  
  display.name    = dataset.phenotype$display.name,
  
  lod.peaks       = list(additive = tibble::as_tibble(dataset.phenotype$lod.peaks$additive),
                         
                         sexInt   = tibble::as_tibble(dataset.phenotype$lod.peaks$sexInt))
  
)

dataset.phenotype <- temp



rm(temp)



library(dplyr)



dataset.phenotype$covar.info <- dataset.phenotype$covar.info[,-2]

dataset.phenotype$covar.info$primary <- FALSE

dataset.phenotype$covar.info$primary[dataset.phenotype$covar.info$sample.column == 'sex'] <- TRUE

dataset.phenotype$covar.info = dataset.phenotype$covar.info[1,]
####



rm(annot_mstrg)
rm(gene_annot_file)
rm(peaks_35peer_qnorm_ALL)
rm(x1,x2,x3,x4)
rm(y1,y2,y3,y4)
rm(Xcovar)
rm(getCisEqtl)
rm(getTransEqtl)
rm(gene_names)
rm(i)
rm(chr)
rm(cis_eqtl_ALL)
rm(trans_eqtl_ALL)
rm(ensembl_names)
rm(annot,annot.mrna,annot.samples,covar.info,covar.matrix,cross_eqtl,data)
rm(lod.peaks)
rm(covar)
rm(peaks_ALLpeer_qnorm_ALL)
rm(perm,perm_a,perm_files,perm_x,pheno_name,pheno_rows,quant_pheno_columns,datatype,display.name)
rm(z,sexInt,raw,qtl_peaks_bin_norm,dataset.protein,dataset.mrna,cross_basic,additive,DO_qtl_int,DO_qtl_int_binary,DO_qtl_scan,DO_qtl_scan_binary,DO_qtl_scan_binary_norm,DO_qtl_scan_normal)
rm(norm_perm,norm_pheno,new_covar,log,k,annot.phenotype,pheno_combined,q,qtl_peaks_binary,qtl_peaks_binary_int,qtl_peaks_both_int,qtl_peaks_both_norm,qtl_peaks_int,qtl_peaks_norm)
save.image("~/Desktop/example/test.RData")

