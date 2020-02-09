
####merge analysis for locus on chromosome 1
###seems like there are two separate associations, TMD, ct.th and ma.ar and ML


set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_REDO.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#load qtl mapping object
load("./results/Rdata/DO_qtl_scan_norm.Rdata")

annot_file = read.csv("./results/flat/annot_file.csv", stringsAsFactors = FALSE)



#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file

##
covar_eqtl = as.matrix(cross_eqtl$covar)
covar_eqtl[,"sex"] = (covar_eqtl[,"sex"] == "M")*1

covar_eqtl = covar_eqtl[,-1]#remove sac date as covar for now

covar_eqtl = apply(covar_eqtl,2,as.numeric) #make sure all cols are numeric
rownames(covar_eqtl) = rownames(cross_eqtl$covar)#make sure rownames match original cross file


##FIX QTL mapping for MAT nonzero vs full
#get qtl list, passed threshold
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

##
####try while Normal####
norm_pheno = as.data.frame(log10(cross_basic$pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,52,54:58,61:70,72,74,76)]))
pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,36,42,50,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}


new_covar = covar
#new_covar[,c(3)] = log10(new_covar[,c(3)])
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.
#####



query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")



#map the phenotypes
#pheno_combined includes normalized and non-normalized phenos, from map_qtl.R
locus1_scan = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)

plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "ML")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar")
plot(locus1_scan, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity")

#we first looked at ML for qsox1. do a snp scan in conf interval
start = 155.05000
end = 155.69592
chr = 1
out_snps_ML_1 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"ML"], k_loco[["1"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                           query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


variants_locus = query_variants(chr, start, end)
genes_locus <- query_genes(chr, start, end)

genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]

if("pseudogene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
}

if("miRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
}

if("rRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
}

top <- top_snps(out_snps_ML_1$lod, out_snps_ML_1$snpinfo)
top[order(top$lod,decreasing = T),]
#top snps are WSB private
plot_snpasso(out_snps_ML_1$lod, out_snps_ML_1$snpinfo, genes = genes_locus)
#condition on top snp rs50769082
snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.4623),
                      sdp=128,
                      snp=c("rs50769082"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs50769082[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs50769082[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus1_scan_cond = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "ML")

plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD")
plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI")
plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax")
plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar")
plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar")
plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar")

plot(locus1_scan_cond, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity")

peaks = find_peaks(locus1_scan_cond, cross_basic$pmap, threshold=4, drop=1.5)

#ML, pMOI,Imax, ct.ar/tt.ar*, tt.ar, ma.ar  and ct.porosity goes away 
#ma.ar drop from 12.8 to 5.9 thresh is 7.6

#TMD doesnt go away but is reduced by about half 

#do snp scan for TMD in CI
start = 154.80370
end = 155.60037
chr = 1

variants_locus = query_variants(chr, start, end)
genes_locus <- query_genes(chr, start, end)

genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]

if("pseudogene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
}

if("miRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
}

if("rRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
}

out_snps_TMD_1 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ct.TMD"], k_loco[["1"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                            query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


top <- top_snps(out_snps_TMD_1$lod, out_snps_TMD_1$snpinfo,drop=20)
top[order(top$lod,decreasing = T),]
plot_snpasso(out_snps_TMD_1$lod, out_snps_TMD_1$snpinfo, genes = genes_locus)

#top snps are also WSB private but seem to be different. kinda two different peaks though. condition on topmost just to see

#condition on top snp rs248974780
snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.0559),
                      sdp=128,
                      snp=c("rs248974780"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}



#redo qtl scan while conditioning on snp
locus1_scan_cond_TMD = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)

plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "ML")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar")
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity")

peaks = find_peaks(locus1_scan_cond_TMD, cross_basic$pmap, threshold=4, drop=1.5)

#everything goes away (except cortical porosity, its just under the lod score (lod 77.37, thresh = 7.77), pos changes a bit (154.817991)! 
#Ok, so if the two snps have the same SDP, why do they have different effects? 
#Maybe this variant is linked to something else?
#look deeper into TMD and ML

out_blup_ML = scan1blup(apr[,1],pheno = pheno_combined[,"ML"], kinship = k_loco[[1]], addcovar =  covar_snp[,c(2,3,4,7:16)],cores = 2)
plot_coefCC(out_blup_ML,map = cross_basic$pmap)
out_blup_TMD = scan1blup(apr[,1],pheno = pheno_combined[,"uCT_Ct.TMD"], kinship = k_loco[[1]], addcovar =  covar_snp[,c(2,3,4,7:16)],cores = 2)
plot_coefCC(out_blup_TMD,map = cross_basic$pmap)

top <- top_snps(out_snps_TMD_1$lod, out_snps_TMD_1$snpinfo,drop = 15)

#rs239189799 = missense variant in Ier5, lod = 30.1
#also multiple missense vars in BC034090, Lod = 28.6

#multiple snps with SDP 192 (WSB and PWK)
#top = 1:155413702_T/C (LOD=22.4)
#condition

snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.413702),
                      sdp=192,
                      snp=c("1:155413702_T/C"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$`A.1:155413702_T/C`[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$`B.1:155413702_T/C`[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus1_scan_cond_155413702 = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)


plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "ML")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar")
plot(locus1_scan_cond_155413702, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity")

#well, they drop, but ML is very near the thresh, and TMD is still there
#maybe because small number of alternate allele?




#what is if condition on both snps that drop TMD but dont ablate it?
snpinfo <- data.frame(chr=c("1","1"),
                      pos=c(155.0559,155.413702),
                      sdp=c(128,192),
                      snp=c("rs248974780","1:155413702_T/C"), stringsAsFactors=FALSE)

snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##

covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 0
  } else{if(covar_snp$B.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 1
  }}
}

covar_snp$alleleB = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$`A.1:155413702_T/C`[i] >0.6){
    covar_snp$alleleB[i] = 1
  } else{if(covar_snp$`B.1:155413702_T/C`[i] >0.6){
    covar_snp$alleleB[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus1_scan_cond_BOTH = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA", "alleleB")],cores = 1)

plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "ML")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar")
plot(locus1_scan_cond_BOTH, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity")

peaks = find_peaks(locus1_scan_cond_BOTH, cross_basic$pmap, threshold=4, drop=1.5)

#they all go away! porosity is suggestive

#do snp scan for TMD in CI
start = 154.803699
end = 155.10491
chr = 1

variants_locus = query_variants(chr, start, end)
genes_locus <- query_genes(chr, start, end)

genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]

if("pseudogene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
}

if("miRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
}

if("rRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
}

out_snps_TMD_bothCOND <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ct.TMD"], k_loco[["1"]],  addcovar =  covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","alleleA","alleleB")],Xcovar=Xcovar,
                            query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)

plot_snpasso(out_snps_TMD_bothCOND$lod, out_snps_TMD_bothCOND$snpinfo, genes = genes_locus)

top <- top_snps(out_snps_TMD_bothCOND$lod, out_snps_TMD_bothCOND$snpinfo,drop = 5)




#compare the two WSB private snps

#condition on top snp rs50769082
snpinfo <- data.frame(chr=c("1","1"),
                      pos=c(155.4623,155.0559),
                      sdp=c(128,128),
                      snp=c("rs50769082","rs248974780"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##

covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs50769082[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs50769082[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

covar_snp$alleleB = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs248974780[i] >0.6){
    covar_snp$alleleB[i] = 1
  } else{if(covar_snp$B.rs248974780[i] >0.6){
    covar_snp$alleleB[i] = 0
  }}
}

which(covar_snp$alleleA != covar_snp$alleleB)
#some difference, some mice are hets for one and full for the other
covar_snp[which(covar_snp$alleleA != covar_snp$alleleB),]

#for example mouse 116
#plot_lodpeaks
DOex = cross_basic
DOex <- DOex[c("1","116"), ]
DOpr <- calc_genoprob(DOex, map=cross_basic$gmap, cores = 1,error_prob = 0.002,map_function="c-f",quiet = F)

DOpr = clean_genoprob(DOpr)

mmarg = maxmarg(DOpr, minprob=0.5)


ph = guess_phase(cross = DOex, geno = mmarg)

plot_onegeno(ph, map = DOex$pmap, ind = "1", chr=1, col=CCcolors)

##
load("~/Desktop/merge_analysis/merge_top_local_eqtl.Rdata")
merge_top_eqtl = merge_top

summary(unlist(lapply(merge_top_eqtl, function(x) nrow(x))))

length(which(unlist(lapply(merge_top_eqtl, function(x) nrow(x))) == 1))

for(i in 1:length(merge_top_eqtl)){
  if("rs248974780" %in% merge_top_eqtl[[i]]$snp_id){
    print(names(merge_top_eqtl)[i])
  }
}
#rs50769082 is an eqtl for Ier5
#so is rs248974780

#plot pxg
plot_lodpeaks(peaks=peaks, map = cross_basic$pmap, chr = 1)
peaks = find_peaks(locus1_scan, cross_basic$pmap, threshold=4, drop=1.5)
#plot lod peaks while scanning for each pheno individually

