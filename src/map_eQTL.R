#map eqtl
library(qtl2)


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

##############
#eQTL mapping

## using only PEER factors (and all of them) (for distal-eqtl) as covars
out_eqtl_distal <- scan1(apr, cross_eqtl$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c(11:ncol(covar))],cores = 20)
save(out_eqtl_distal,file ="out_eqtl_distal.Rdata")

distal_eqtl_peaks = find_peaks(out_eqtl_distal, cross_basic$pmap, threshold=4, drop=1.5)

#using first 35 PEER factors and sex as covars (local eQTL)
out_eqtl_local <- scan1(apr, cross_eqtl$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,11:45)],cores = 20)
save(out_eqtl_local,file ="out_eqtl_local.Rdata")

local_qtl_peaks = find_peaks(out_eqtl_local, cross_basic$pmap, threshold=4, drop=1.5)

