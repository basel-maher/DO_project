#QTL mapping 
library(qtl2)

#load the allele probs
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
#load(file = "./results/Rdata/k_basic_cleaned.Rdata") #overall (not used)
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO

#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file

#qtl mapping for continuous traits
#which traits require FL as covar? bending traits?
DO_qtl_scan = scan1(apr, cross_basic$pheno[,c(5:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])
save(DO_qtl_scan,file = "DO_qtl_scan.Rdata")
#Warning: Can't fit binary model with kinship matrix
DO_qtl_scan_binary = scan1(apr, cross_combined$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])
save(DO_qtl_scan_binary,file = "DO_qtl_scan_binary.Rdata")

qtl_peaks = find_peaks(DO_qtl_scan, cross_combined$pmap, threshold=4, drop=1.5)
qtl_peaks_binary = find_peaks(DO_qtl_scan_binary, cross_combined$pmap, threshold=4, drop=1.5)

qtl_peaks_both = rbind(qtl_peaks,qtl_peaks_binary)
write.csv(qtl_peaks_both, file = "qtl_peaks.csv",row.names = FALSE,quote = FALSE)
