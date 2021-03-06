#QTL mapping 
set.seed(8675309)
library(qtl2)
#load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
#load the allele probs
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
#load(file = "./results/Rdata/k_basic_cleaned.Rdata") #overall (not used)
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file


###

#normalize  phenos##

quant_pheno_columns = c(5:70,72,74,76)#quantitative phenotypes

for(i in quant_pheno_columns){
  print(paste(i,shapiro.test(cross_basic$pheno[,i])$p.value)) #if pval < alpha, not Normal 
  
  if(shapiro.test(cross_basic$pheno[,i])$p.value>0.05){
    print(i)
  }
}
#only 15,18,19,20,22,36,42,50,53,59, 60 are Normal



#qtl mapping for continuous traits

DO_qtl_scan = scan1(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
#save(DO_qtl_scan,file = "./results/Rdata/DO_qtl_scan.Rdata")

#scan MAT as binary traits
DO_qtl_scan_binary = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
#save(DO_qtl_scan_binary,file = "./results/Rdata/DO_qtl_scan_binary.Rdata")
load("./results/Rdata/DO_qtl_scan_binary.Rdata")

#find peaks and bind them together
qtl_peaks = find_peaks(DO_qtl_scan, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_binary = find_peaks(DO_qtl_scan_binary, cross_basic$pmap, threshold=4, drop=1.5)

qtl_peaks_both = rbind(qtl_peaks,qtl_peaks_binary)
write.csv(qtl_peaks_both, file = "./results/flat/qtl_peaks.csv",row.names = FALSE,quote = FALSE)


####try after transforming####
norm_pheno = as.data.frame(cross_basic$pheno)

norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1 + 1
norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2 + 1
norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1

norm_pheno$bending_work_post_yield = norm_pheno$bending_work_post_yield + 1
norm_pheno$bending_PYD = norm_pheno$bending_PYD + 1

norm_pheno = as.data.frame(log10(norm_pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,52,54:58,61:70,72,74,76)]))

pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,36,42,50,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.


DO_qtl_scan_normal = scan1(apr, pheno_combined, k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
#save(DO_qtl_scan_normal,file = "./results/Rdata/DO_qtl_scan_norm.Rdata")
load("./results/Rdata/DO_qtl_scan_norm.Rdata")

qtl_peaks_norm = find_peaks(DO_qtl_scan_normal, cross_basic$pmap, threshold=4, drop=1.5)
##cant account for kinship in binary model
DO_qtl_scan_binary_norm = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
#save(DO_qtl_scan_binary_norm,file = "./results/Rdata/DO_qtl_scan_bin_norm.Rdata")
load("./results/Rdata/DO_qtl_scan_bin_norm.Rdata")


qtl_peaks_bin_norm = find_peaks(DO_qtl_scan_binary_norm, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_both_norm = rbind(qtl_peaks_norm,qtl_peaks_bin_norm)


#
#
#

#calc heritability
h = est_herit(pheno = pheno_combined, kinship = k,addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")] )

h_df = as.data.frame(h)
#
#
#


# 
# ### interactive model w/ sex
# DO_qtl_int = scan1(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = covar[,"sex"],cores = 2)
# #save(DO_qtl_int,file = "./results/Rdata/DO_qtl_int.Rdata")
# load("./results/Rdata/DO_qtl_int.Rdata")
# qtl_peaks_int = find_peaks(DO_qtl_int, cross_basic$pmap, threshold=4, drop=1.5)
# 
# 
# DO_qtl_int_binary = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = covar[,"sex"],cores = 2)
# #save(DO_qtl_int_binary,file = "./results/Rdata/DO_qtl_int_binary.Rdata")
# load("./results/Rdata/DO_qtl_int_binary.Rdata")
# 
# qtl_peaks_binary_int = find_peaks(DO_qtl_int_binary, cross_basic$pmap, threshold=4, drop=1.5)
# qtl_peaks_both_int = rbind(qtl_peaks_int,qtl_peaks_binary_int)
# 
# 
# ##diff both
# x = merge(qtl_peaks_both_int,qtl_peaks_both_norm, by=c("lodcolumn","chr"))
# 
# #same but norm
# DO_qtl_int_norm = scan1(apr, pheno_combined, k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = new_covar[,"sex"],cores = 2)
# save(DO_qtl_int_norm,file = "./results/Rdata/DO_qtl_int_norm.Rdata")
# load("./results/Rdata/DO_qtl_int_norm.Rdata")
# 
# qtl_peaks_int_norm = find_peaks(DO_qtl_int_norm, cross_basic$pmap, threshold=4, drop=1.5)
# 
# DO_qtl_int_binary_norm = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = new_covar[,"sex"],cores = 2)
# save(DO_qtl_int_binary_norm,file = "./results/Rdata/DO_qtl_int_binary_norm.Rdata")
# load("./results/Rdata/DO_qtl_int_binary_norm.Rdata")
# qtl_peaks_int_bin_norm = find_peaks(DO_qtl_int_binary_norm, cross_basic$pmap, threshold=4, drop=1.5)
# 
# qtl_peaks_both_int = rbind(qtl_peaks_int_bin_norm,qtl_peaks_int_norm)
# 
# 
# 
# 


#use qtl2::scan1perm() for permutations. We did 1000 permutations for each pheno.
#this was done on our supercomputing cluster
#example:
#scan1perm(apr,pheno_combined,k_loco, Xcovar=Xcovar, addcovar = cov[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],n_perm = 1000,perm_Xsp = TRUE,chr_lengths=chr_lengths(cross_basic$gmap)
#

#####get qtl that pass perm threshold
qtl_peaks_both$perm_thresh = NA

perm = list.files("./results/Rdata/qtl_perms/")
perm_files = perm[grep("^perms_",perm)]

for(i in 1:length(perm_files)){
  
  load(paste0("./results/Rdata/qtl_perms/",perm_files[i]))
  perm_a = summary(perm)$A[1]
  perm_x = summary(perm)$X[1]
  
  pheno_name = perm_files[i]
  pheno_name = gsub(x = perm_files[i],pattern = "perms_",replacement = "")
  pheno_name = gsub(x = pheno_name,pattern = ".Rdata",replacement = "")
  
  pheno_rows = which(qtl_peaks_both$lodcolumn == pheno_name)
  
  for(i in 1:length(pheno_rows)){
    if(qtl_peaks_both$chr[pheno_rows[i]] == "X"){
      qtl_peaks_both$perm_thresh[[pheno_rows[i]]] = perm_x
    } else {qtl_peaks_both$perm_thresh[[pheno_rows[i]]] = perm_a}
  }
}

qtl = qtl_peaks_both[which(qtl_peaks_both$lod >= qtl_peaks_both$perm_thresh),]




##same for norm
qtl_peaks_both_norm$perm_thresh = NA
perm = list.files("./results/Rdata/qtl_perms/")

perm_files = perm[grep("norm_perms_",perm)]
#remove INT
perm_files = perm_files[-grep("_INT_",perm_files)]

for(i in 1:length(perm_files)){
  
  load(paste0("./results/Rdata/qtl_perms/",perm_files[i]))

  perm_a = summary(norm_perm)$A[1]
  perm_x = summary(norm_perm)$X[1]
  
  pheno_name = perm_files[i]
  pheno_name = gsub(x = perm_files[i],pattern = "norm_perms_",replacement = "")
  pheno_name = gsub(x = pheno_name,pattern = ".Rdata",replacement = "")
  
  pheno_rows = which(qtl_peaks_both_norm$lodcolumn == pheno_name)
  
  for(i in 1:length(pheno_rows)){
    if(qtl_peaks_both_norm$chr[pheno_rows[i]] == "X"){
      qtl_peaks_both_norm$perm_thresh[[pheno_rows[i]]] = perm_x
    } else {qtl_peaks_both_norm$perm_thresh[[pheno_rows[i]]] = perm_a}
  }
}
qtl_norm = qtl_peaks_both_norm[which(qtl_peaks_both_norm$lod >= qtl_peaks_both_norm$perm_thresh),]
write.csv(qtl_norm, file = "./results/flat/qtl_norm_pass_thresh", row.names = FALSE, quote = FALSE)
###


# #same for NORM INT
# qtl_peaks_both_int$perm_thresh = NA
# perm = list.files("./results/Rdata/qtl_perms/")
# 
# perm_files = perm[grep("norm_perms_INT_",perm)]
# 
# for(i in 1:length(perm_files)){
#   
#   load(paste0("./results/Rdata/qtl_perms/",perm_files[i]))
#   
#   perm_a = summary(norm_perm)$A[1]
#   perm_x = summary(norm_perm)$X[1]
#   
#   pheno_name = perm_files[i]
#   pheno_name = gsub(x = perm_files[i],pattern = "norm_perms_INT_",replacement = "")
#   pheno_name = gsub(x = pheno_name,pattern = ".Rdata",replacement = "")
#   
#   pheno_rows = which(qtl_peaks_both_int$lodcolumn == pheno_name)
#   
#   for(i in 1:length(pheno_rows)){
#     if(qtl_peaks_both_int$chr[pheno_rows[i]] == "X"){
#       qtl_peaks_both_int$perm_thresh[[pheno_rows[i]]] = perm_x
#     } else {qtl_peaks_both_int$perm_thresh[[pheno_rows[i]]] = perm_a}
#   }
# }
# 
# qtl_norm_INT = qtl_peaks_both_int[which(qtl_peaks_both_int$lod >= qtl_peaks_both_int$perm_thresh),]
# write.csv(qtl_norm_INT, file = "./results/flat/qtl_norm_INT_pass_thresh", row.names = FALSE, quote = FALSE)
# 
# ###


#full_qtl = merge(qtl_norm, qtl_norm_INT, by=c("lodcolumn","chr"))
#full_qtl$LODi = full_qtl$lod.y - full_qtl$lod.x

# NOT USED
# PLOTTING, ETC
# 
# 
# 
# TMD_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Ct.TMD"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(TMD_blup, file = "./results/Rdata/TMD_blup.Rdata")
# 
# #load("./results/Rdata/TMD_blup.Rdata")
# png("~/Desktop/TMD.png")
# plot_coefCC(TMD_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=45), main = "TMD")
# dev.off() 
# 
# 
# plot_coefCC(TMD_blup, cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=45), main = "TMD")
# 
# TMD_coef_scan = scan1coef(apr,pheno = pheno_combined[,"uCT_Ct.TMD"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# 
# pMOI_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_pMOI"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(pMOI_blup, file = "./results/Rdata/pMOI_blup.Rdata")
# load("./results/Rdata/pMOI_blup.Rdata")
# 
# png("~/Desktop/pMOI.png")
# plot_coefCC(pMOI_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=47), main = "pMOI")
# dev.off()
# 
# imax_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Imax"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(imax_blup, file = "./results/Rdata/imax_blup.Rdata")
# load("./results/Rdata/imax_blup.Rdata")
# 
# png("~/Desktop/imax.png")
# plot_coefCC(imax_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=48), main = "Imax")
# dev.off()
# 
# ttAr_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Tt.Ar"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(ttAr_blup, file = "./results/Rdata/ttAr_blup.Rdata")
# load("./results/Rdata/ttAr_blup.Rdata")
# 
# png("~/Desktop/ttar.png")
# plot_coefCC(ttAr_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=43), main = "Tt.Ar")
# dev.off()
# 
# maAr_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Ma.Ar"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(maAr_blup, file = "./results/Rdata/maAr_blup.Rdata")
# load("./results/Rdata/maAr_blup.Rdata")
# 
# png("~/Desktop/Ma.Ar.png")
# plot_coefCC(maAr_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=42), main = "Ma.Ar")
# dev.off()
# 
# ML_blup = scan1blup(apr,pheno = pheno_combined[,"ML" ], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(ML_blup, file = "./results/Rdata/ML_blup.Rdata")
# load("./results/Rdata/ML_blup.Rdata")
# 
# png("~/Desktop/ml.png")
# plot_coefCC(ML_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=10),main = "ML")
# dev.off()
# 
# porosity_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Ct.porosity"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# save(porosity_blup, file = "./results/Rdata/porosity_blup.Rdata")
# load("./results/Rdata/porosity_blup.Rdata")
# 
# png("~/Desktop/porosity.png")
# plot_coefCC(porosity_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=46),main="Ct.Porosity")
# dev.off()
# 
# #
# TMD_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
# TMD_coef = TMD_blup[TMD_mark,1:8]
# 
# pMOI_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
# pMOI_coef = pMOI_blup[pMOI_mark,1:8]
# 
# imax_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
# imax_coef = imax_blup[imax_mark,1:8]
# 
# ttAr_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.19293)
# ttAr_coef = ttAr_blup[ttAr_mark,1:8]
# 
# maAr_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.32971)
# maAr_coef = maAr_blup[maAr_mark,1:8]
# 
# 
# ML_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.35757)
# ML_coef = ML_blup[ML_mark,1:8]
# 
# por_mark = find_marker(cross_basic$pmap, chr = 1, pos = 156.18734)
# por_coef = porosity_blup[por_mark,1:8]
# 
# cor(pMOI_coef, TMD_coef, method = "k")




###
#make qtl_loci
#qtl file
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

#remove FFP and soleus
qtl_norm = qtl_norm[-c(1:3),]
#remove MAT_vol1_nonzero
qtl_norm = qtl_norm[-(which(qtl_norm$lodcolumn == "MAT_VOL1_nonzero")),]

#define loci
qtl_norm$locus = 0

loc_idx = 1
qtl_loc = as.data.frame(matrix(nrow=nrow(qtl_norm), ncol=ncol(qtl_norm)))
colnames(qtl_loc) = colnames(qtl_norm)
for (i in c(1,2,3,4,8,10,16,"X")){
  sub = subset(qtl_norm, qtl_norm$chr == i)
  
  while(any(sub$locus == 0)){
    min_sub = min(sub$pos)
    idx = which((sub$pos >= min_sub) & (sub$pos <= min_sub + 1.5))
    sub[idx,"locus"] = loc_idx
    qtl_loc = rbind(qtl_loc, sub[idx,])
    sub = sub[-idx,]
    loc_idx = loc_idx + 1
  }
}

qtl_loc = qtl_loc[-which(is.na(qtl_loc)),]

write.csv(qtl_loc, file = "./results/flat/qtl_loc", quote = FALSE,row.names = FALSE)







