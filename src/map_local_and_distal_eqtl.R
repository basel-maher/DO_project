#Map distal and local eQTL

library(qtl2)
library(GenomicRanges)
library(rtracklayer)

#library(biomaRt)
#for mapping, you need the allele probs, kinship and cross file
load("./results/Rdata/apr_basic_cleaned.Rdata")
load("./results/Rdata/cross_eqtl_REDO.Rdata")
load("./results/Rdata/k_loco_basic_cleaned.Rdata")



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
###########
#This function allows you to merge a peak file with the 
#output of the biomaRt function below, which would give you the TSS positions and IDs for all genes in the peakfile.
#This function then return a dataframe with trans_eqtls given the arguments. In this case, trans eQTL are
#defined as being on a different chrom OR same chrom but outside the given transMinDistFromTSS (greater than or equal to this minimum distance)

getDistalEqtl = function(peaks,geneAnnot,lodThreshAuto,lodThreshX, distalMinDistFromStart,geneCol1,geneCol2){
  out = merge(peaks,geneAnnot,by.x = geneCol1, by.y = geneCol2)
  
  out$dist_start = abs((out$pos*1000000) - out$Start) #multiply by a mil to get bp, calculate distance from start
  
  idx = which((out$chr != "X") & (out$lod >= lodThreshAuto)) #find autosomal peaks that meet threshold
  idx = c(idx, which((out$chr == "X") & (out$lod >= lodThreshX))) #find X peaks that meet threshold
  out = out[idx,] #must meet or exceed lod thresh
  
  out$distal = 0 #set col of zeros
  out$distal[which(as.numeric(out$dist_start) >=distalMinDistFromStart)] = 1 #if meets distance criterion for distal eqtl, make = 1
  
  out$dist_start[which(out$chr != out$Reference)] = "diff_chrom" #if not on same chromosome, make = "diff_chrom"
  out$distal[which(out$dist_start=="diff_chrom")] = 1 #make = 1 if on a different chrom
  #
  out = out[which(out$distal == 1),] #subset to only "distal" eqtl
  #out = out[,c(1:11)] # take these columns
  out = out[which(out$chr %in% c(seq(1,19,1),"X")),] #take those only on chroms 1:19 and X 
  return(out)
}
#########
#########
#########
#############################################################################################################################

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

#using first sex and all 48 PEER factors and as covars (local eQTL)
out_eqtl_local <- scan1(apr, cross_eqtl$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,11:58)],cores = 20)

local_eqtl_peaks = find_peaks(out_eqtl_local, cross_eqtl$pmap, threshold=4, drop=1.5)

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
annot_file = read.csv("./results/flat/annot_file.csv", stringsAsFactors = FALSE)

#annot_file = readGFF("./results/flat/RNA-seq/mus_stringtie_merged.gtf")
#annot_file = annot_file[which(annot_file$type=="transcript"),]
#remove transcripts not on somatic and X chroms
chr = c(seq(1:19),"X")
annot_file = annot_file[which(annot_file$Reference %in% chr),]


#annot_file = annot_file[,c(9,11)]
#annot_file = unique(annot_file)

#CHANGE TO HAVE 1 ENTRY PER MSTRG/ID
#temp
#annot_file = annot_file[-which(duplicated(annot_file$gene_id)),]
# 

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

distal_eqtl = getDistalEqtl(peaks,annot_file,lodThreshAuto = 10.89, lodThreshX = 11.55, distalMinDistFromStart = 1000000, geneCol1 = "lodcolumn", geneCol2 = "Gene.ID")

###remove -ps pseudogenes, -rs and GM genes from trans eqtl list
##(what about -rs? related sequence)
distal_eqtl = distal_eqtl[-grep("-ps",distal_eqtl$Gene.Name),]
distal_eqtl = distal_eqtl[-grep("-rs",distal_eqtl$Gene.Name),]
distal_eqtl = distal_eqtl[-grep("Gm[0-9]",distal_eqtl$Gene.Name),]

save(distal_eqtl, file = "./results/Rdata/distal_eqtl.Rdata")
##CHANGE TO CALCULATE USING MSTRG INSTEAD OF SYMBOL





#####
#calc distal eqtl hotspots
######trans_eQTL hotspots######
#convert the genome into intervals
#get ensembl assembly index .fai file
destfile = "~/Desktop/DO_proj/data/Mus_musculus.GRCm38.dna.toplevel.fa.gz.fai"
download.file("ftp://ftp.ensembl.org/pub/release-95/fasta/mus_musculus/dna_index/Mus_musculus.GRCm38.dna.toplevel.fa.gz.fai",destfile = destfile)

mus_fai = read.delim(destfile,header = FALSE,col.names = c("chrom","length","x1","x2","x3"))

hotspot_intervals = as.data.frame(matrix(nrow = 1,ncol = 3))
colnames(hotspot_intervals) = c("int_start","int_end","chr")

for(i in 1:21){
  mod = mus_fai$length[i] %% 4000000 
  df = as.data.frame(matrix(seq(1, mus_fai$length[i],by=4000000)))
  df$interval_end = df$V1+4000000-1
  
  if(mod != 0){
    df[dim(df)[1],2] = mus_fai$length[i]
  }
  
  df$chr = mus_fai[i,1]
  
  colnames(df) = colnames(hotspot_intervals)
  hotspot_intervals = rbind.data.frame(hotspot_intervals,df)
  
  df = as.data.frame(matrix(seq(1000000, mus_fai$length[i],by=4000000)))
  df$interval_end = df$V1+4000000-1
  
  if(mod != 0){
    df[dim(df)[1],2] = mus_fai$length[i]
  }
  
  df$chr = mus_fai[i,1]
  
  colnames(df) = colnames(hotspot_intervals)
  hotspot_intervals = rbind.data.frame(hotspot_intervals,df)
  
  df = as.data.frame(matrix(seq(2000000, mus_fai$length[i],by=4000000)))
  df$interval_end = df$V1+4000000-1
  
  if(mod != 0){
    df[dim(df)[1],2] = mus_fai$length[i]
  }
  
  df$chr = mus_fai[i,1]
  
  colnames(df) = colnames(hotspot_intervals)
  hotspot_intervals = rbind.data.frame(hotspot_intervals,df)
  
  df = as.data.frame(matrix(seq(3000000, mus_fai$length[i],by=4000000)))
  df$interval_end = df$V1+4000000-1
  
  if(mod != 0){
    df[dim(df)[1],2] = mus_fai$length[i]
  }
  
  df$chr = mus_fai[i,1]
  
  colnames(df) = colnames(hotspot_intervals)
  hotspot_intervals = rbind.data.frame(hotspot_intervals,df)
  
  
}




hotspot_intervals = rbind(hotspot_intervals,c(1,mus_fai$length[22],"MT"))#can add MT
hotspot_intervals = hotspot_intervals[-1,]
hotspot_intervals$n_trans = 0


distal_eqtl$pos = as.numeric(distal_eqtl$pos*1000000)
###
for(i in 1:nrow(hotspot_intervals)){
  print(i)
  chr = hotspot_intervals[i,3]
  chr = sub(pattern="chr", x=chr, replacement="")
  sub = distal_eqtl[which(distal_eqtl$chr == chr),]
  
  start_int = as.numeric(hotspot_intervals[i,1])
  end_int = as.numeric(hotspot_intervals[i,2])
  
  num_trans = length(which(sub$pos >=start_int & sub$pos<= end_int))
  hotspot_intervals$n_trans[i] = num_trans
  #summary(distal_eqtl$ci_hi - distal_eqtl$ci_lo) remove stupid large CIs
}

#hotspot_intervals$chr = paste0("chr",hotspot_intervals$chr)

###hotspot interals permutation. For each interval, count number of 
##copy trans_eqtl_ALL, make chr and pos random
vec = c()
chr_vec = c(1:19,"X")
hot = hotspot_intervals

for(i in 1:1000){
  print(i)
  tab = distal_eqtl[,c(1,2,3)]
  tab[,c(2,3)]=NA
  #hot = hotspot_intervals
  int = i+4
  hot[,int] = NA
  for(j in 1:nrow(tab)){
    tab[j,2] = sample(chr_vec,1)
    tab[j,3] = sample(1:mus_fai[which(mus_fai$chrom == sub(pattern="chr", x=tab[j,2], replacement="")),"length"],1)
  }
  #tab$chr = paste0("chr",tab$chr)
  for(k in 1:nrow(hot)){
    #print(k)
    chr = hotspot_intervals[k,3]
    sub = tab[which(tab$chr == chr),]
    
    start_int = as.numeric(hot[k,1])
    end_int = as.numeric(hot[k,2])
    
    num_trans = length(which(sub$pos >=start_int & sub$pos<= end_int))
    hot[k,int] = num_trans
    #summary(distal_eqtl$ci_hi - distal_eqtl$ci_lo) remove stupid large CIs
  }
  #vec = c(vec,hot$n_trans)
}

hot$mean = apply(hot,1,function(x) mean(as.numeric(x[5:1004])))
hot$median = apply(hot,1,function(x) median(as.numeric(x[5:1004])))
hot$sd = apply(hot,1,function(x) sd(as.numeric(x[5:1004])))
hot$max = apply(hot,1,function(x) max(as.numeric(x[5:1004])))
hot$min = apply(hot,1,function(x) min(as.numeric(x[5:1004])))

hot = hot[,-c(5:1004)]











