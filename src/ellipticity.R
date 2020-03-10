set.seed(8675309)
library(qtl2)


#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO

#get Xcovar
Xcovar <- get_x_covar(cross_basic)


###
#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file



#calculate ellipticity
a = cross_basic$pheno[,"uCT_Imax"]
b = cross_basic$pheno[,"uCT_Imin"]

e = sqrt((a^2 - b^2)/a^2)


#map ellipticity
ellipticity_scan = scan1(apr, e, k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])
e_peaks = find_peaks(ellipticity_scan, cross_basic$pmap, threshold=4, drop=1.5)

