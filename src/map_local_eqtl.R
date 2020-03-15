#Map distal and local eQTL

library(qtl2)
library(GenomicRanges)
library(rtracklayer)

#library(biomaRt)
#for mapping, you need the allele probs, kinship and cross file
load("./results/Rdata/apr_basic_cleaned.Rdata")
load("./results/Rdata/cross_eqtl_REDO.Rdata")
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



#########################################################################################################################################
#                                                         FUNCTIONS
#########################################################################################################################################
#define functions for calculating distal and local eqtl
#These use lod thresholds and a distance from the TSS to define distal and local eQTL peaks
#return dataframes with local or distal eqtl

getLocalEqtl = function(peaks,geneAnnot,lodThreshAuto, lodThreshX, localDistFromStart,geneCol1,geneCol2){
  out = merge(peaks,geneAnnot,by.x = geneCol1, by.y = geneCol2)
  
  out$dist_start = abs((out$pos*1000000) - out$Start) #multiply by a mil to get bp, calculate distance from start
  out = out[which(out$chr == out$Reference),] #must be on same chrom
  out = out[which(out$dist_start <=localDistFromStart),] #must be within a certain distance from start
  idx = which((out$chr != "X") & (out$lod >= lodThreshAuto)) #find autosomal peaks that meet threshold
  idx = c(idx, which((out$chr == "X") & (out$lod >= lodThreshX))) #find X peaks that meet threshold
  out = out[idx,] #must meet or exceed lod thresh
  return(out)
}


##############
#eQTL mapping
#This was done on our high performance computing cluster. 
#We split our transcripts into 5 pools and ran each independently, in order to speed up calculations.

#using first sex and all 48 PEER factors and as covars (local eQTL)
#out_eqtl_local <- scan1(apr, cross_eqtl$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,11:58)],cores = 20)

#local_eqtl_peaks = find_peaks(out_eqtl_local, cross_eqtl$pmap, threshold=4, drop=1.5)



annot_file = read.csv("./results/flat/annot_file.csv", stringsAsFactors = FALSE)


chr = c(seq(1:19),"X")
annot_file = annot_file[which(annot_file$Reference %in% chr),]



p1 = read.csv("./results/flat/RNA-seq/eqtl_peaks_REDO/local_eqtl_peaks_1.csv",stringsAsFactors = FALSE)
p2 = read.csv("./results/flat/RNA-seq/eqtl_peaks_REDO/local_eqtl_peaks_2.csv",stringsAsFactors = FALSE)
p3 = read.csv("./results/flat/RNA-seq/eqtl_peaks_REDO/local_eqtl_peaks_3.csv",stringsAsFactors = FALSE)
p4 = read.csv("./results/flat/RNA-seq/eqtl_peaks_REDO/local_eqtl_peaks_4.csv",stringsAsFactors = FALSE)
p5 = read.csv("./results/flat/RNA-seq/eqtl_peaks_REDO/local_eqtl_peaks_5.csv",stringsAsFactors = FALSE)
peaks = rbind(p1,p2,p3,p4,p5)
peaks = peaks[,-1]

#prune output to only include only those that pass LOD threshold from "./src/calc_eqtl_perms.R"

local_eqtl = getLocalEqtl(peaks,annot_file,lodThreshAuto = 10.89, lodThreshX = 11.55, localDistFromStart = 1000000, geneCol1 = "lodcolumn", geneCol2 = "Gene.ID")
save(local_eqtl, file = "./results/Rdata/local_eqtl.Rdata")
#

