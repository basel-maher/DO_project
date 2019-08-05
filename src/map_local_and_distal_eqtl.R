#Map distal and local eQTL

library(qtl2)
library(biomaRt)
#for mapping, you need the allele probs, kinship and cross file
load("./results/Rdata/apr_basic_cleaned.Rdata")
load("./results/Rdata/paper_analyses/cross_eqtl.Rdata")
load("./results/Rdata/k_loco_basic_cleaned.Rdata")

#get the X chrom covars from the cross file
Xcovar <- get_x_covar(cross_eqtl)


#create a covar object from covariates in cross file
#must be numeric
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1 #convert sex to 1's and 0's
covar[,1] = as.factor(covar[,1]) #sac date to factors
covar[,6] = as.factor(covar[,6]) #generation to factors

covar = apply(covar,2,as.numeric)
rownames(covar) = rownames(cross_eqtl$covar)

##############
#eQTL mapping
#This was done on our high performance computing cluster. 
#We split our transcripts into 5 pools and ran each independently, in order to speed up calculations.

#using first 35 PEER factors and sex as covars (local eQTL)
out_eqtl_local <- scan1(apr, cross_eqtl$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,17:51)],cores = 20)

local_eqtl_peaks = find_peaks(out_eqtl_local, cross_eqtl$pmap, threshold=4, drop=1.5)
local_eqtl_peaks = local_eqtl_peaks[,-2]

#convert gene names to symbol using transcript quant output
annot = read.delim("./data/314-FarberDO2_S6.gene_abund.tab", stringsAsFactors = FALSE)#read annotation file

local_eqtl_peaks$lodcolumn = annot[match(local_eqtl_peaks$lodcolumn,annot$Gene.ID),"Gene.Name"] #match ID to get gene name
write.csv(local_eqtl_peaks, file = "./results/flat/local_eqtl_peaks.csv", quote=FALSE)

#using 48 PEER factors as covariates (distal eQTL)
out_eqtl_distal <- scan1(apr, cross_eqtl$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c(17:ncol(covar))],cores = 20)

distal_eqtl_peaks = find_peaks(out_eqtl_distal, cross_eqtl$pmap, threshold=4, drop=1.5)
distal_eqtl_peaks = distal_eqtl_peaks[,-2]

#convert gene names to symbol using transcript quant output
annot = read.delim("./data/314-FarberDO2_S6.gene_abund.tab", stringsAsFactors = FALSE)#read annotation file

distal_eqtl_peaks$lodcolumn = annot[match(distal_eqtl_peaks$lodcolumn,annot$Gene.ID),"Gene.Name"] #match ID to get gene name

write.csv(distal_eqtl_peaks, file = "./results/flat/distal_eqtl_peaks.csv", quote=FALSE)

###
#generate gene annotation file. This contains ensembl gene IDs, chromosomes, and TSS start and end sites for all genes in our peak data
#uses biomaRt
#listMarts() indicates ENSEMBL version 97 is used = GRCm38 (mm10)
#GigaMUGA data is GRCm38. Same as HISAT2 genome_snp

# genes = unique(c(local_eqtl_peaks$lodcolumn, distal_eqtl_peaks$lodcolumn))
# 
# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# 
# biomart_genes = getBM(attributes = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position",'transcript_end'), mart = mouse)[,c(1:5)]#external_gene_name)
# 
# biomart_genes = biomart_genes[which(biomart_genes$mgi_symbol %in% genes),] #prune to contain only genes in our data
# 
# #which of our genes cannot be found in biomart?
# no_mart = genes[which(genes %in% biomart_genes$mgi_symbol == FALSE)]
# #annotation file from quantification
# annot = read.delim("~/Desktop/DO_proj/data/networks/314-FarberDO2_S6.gene_abund.tab")
# 
# #use quantification file to annotate genes that cant be found in biomart
# gene_annot_file = annot[which(annot$Gene.Name %in% no_mart),]
# gene_annot_file = gene_annot_file[,c(1,2,3,5,6)]
# 
# colnames(gene_annot_file) = colnames(biomart_genes)
# gene_annot_file = rbind(gene_annot_file,biomart_genes)#bind the annotations from biomart with the non-biomart ones
# 
# gene_annot_file = gene_annot_file[-which(duplicated(gene_annot_file$mgi_symbol)),]

#USE ONLY ANNOTATION FILE FROM RNASEQ QUANTIFICATION
genes = unique(c(local_eqtl_peaks$lodcolumn, distal_eqtl_peaks$lodcolumn))

gene_annot_file = annot[which(annot$Gene.Name %in% genes),] #use only mapped genes
gene_annot_file = gene_annot_file[,c(1,2,3,5,6)] #take useful columns

colnames(gene_annot_file) = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position" ,"transcript_end")
write.csv(gene_annot_file, file = "./results/flat/gene_annot_file.csv", quote = FALSE, row.names = FALSE)
# 





#prune output to only include only those that pass LOD threshold from "./src/calc_eqtl_perms.R"

#define functions for calculating distal and local eqtl
#These use lod thresholds and a distance from the TSS to define distal and local eQTL peaks
#return dataframes with local or distal eqtl

getLocalEqtl = function(peaks,geneAnnot,lodThresh,localDistFromStart,geneCol1,geneCol2){
  out = merge(peaks,geneAnnot,by.x = geneCol1, by.y = geneCol2)

  out$dist_start = abs((out$pos*1000000) - out$start_position) #multiply by a mil to get bp, calculate distance from start
  out = out[which(out$chr == out$chromosome_name),] #must be on same chrom
  out = out[which(out$dist_start <=localDistFromStart),] #must be within a certain distance from start
  out = out[which(out$lod>=lodThresh),] #must meet or exceed lod thresh
  return(out)
}
###########
#This function allows you to merge a peak file with the 
#output of the biomaRt function below, which would give you the TSS positions and IDs for all genes in the peakfile.
#This function then return a dataframe with trans_eqtls given the arguments. In this case, trans eQTL are
#defined as being on a different chrom OR same chrom but outside the given transMinDistFromTSS (greater than or equal to this minimum distance)

getDistalEqtl = function(peaks,geneAnnot,lodThresh,distalMinDistFromStart,geneCol1,geneCol2){
  out = merge(peaks,geneAnnot,by.x = geneCol1, by.y = geneCol2)

  out$dist_start = abs((out$pos*1000000) - out$start_position) #multiply by a mil to get bp, calculate distance from start
  out = out[which(out$lod>=lodThresh),]
  out$distal = 0 #set col of zeros
  out$distal[which(out$chr != out$chromosome_name)] = 1 #if not on same chromosome, make = 1
  out$dist_start[which(out$distal==1)] = "diff_chrom" #annotate as being on a different chromosome
  #
  out$distal[which(as.numeric(out$dist_start) >=distalMinDistFromStart)] = 1 #if meets distance criterion for distal eqtl, make = 1
  out = out[which(out$distal == 1),] #subset to only "distal" eqtl
  out = out[,c(1:13)] # take these columns
  out = out[which(out$chromosome_name %in% c(seq(1,19,1),"X")),] #take those only on chroms 1:19 and X
  return(out)
}







local_autosomal = 9.9
local_X = 10.9

local_idx = which((local_eqtl_peaks$chr != "X") & (local_eqtl_peaks$lod >= local_autosomal))
local_idx = c(local_idx, which((local_eqtl_peaks$chr == "X") & (local_eqtl_peaks$lod >= local_X)))

x = 
distal_autosomal = 10.7
distal_X = 11

distal_idx = which((distal_eqtl_peaks$chr != "X") & (distal_eqtl_peaks$lod >= distal_autosomal))
distal_idx = c(distal_idx, which((distal_eqtl_peaks$chr == "X") & (distal_eqtl_peaks$lod >= distal_X)))
