#QTL mapping 
library(qtl2)
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
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
}
#only 15,18,19,20,22,33,36,42,50,53,59 are Normal



#qtl mapping for continuous traits

#which traits require FL as covar? bending traits?
DO_qtl_scan = scan1(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(DO_qtl_scan,file = "./results/Rdata/DO_qtl_scan.Rdata")


DO_qtl_scan_binary = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(DO_qtl_scan_binary,file = "./results/Rdata/DO_qtl_scan_binary.Rdata")

#find peaks and bind them together
qtl_peaks = find_peaks(DO_qtl_scan, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_binary = find_peaks(DO_qtl_scan_binary, cross_basic$pmap, threshold=4, drop=1.5)

qtl_peaks_both = rbind(qtl_peaks,qtl_peaks_binary)
write.csv(qtl_peaks_both, file = "qtl_peaks.csv",row.names = FALSE,quote = FALSE)


####try while Normal####
norm_pheno = as.data.frame(log10(cross_basic$pheno[,c(6:14,16,17,21,23:32,34,35,37:41,43:49,51,52,54:58,60:70,72,74,76)])+0.0000001)
pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,33,36,42,50,53,59)])
new_covar = covar
new_covar[,c(3)] = log10(new_covar[,c(3)])+0.0000001

DO_qtl_scan_normal = scan1(apr, pheno_combined, k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(DO_qtl_scan_normal,file = "./results/Rdata/DO_qtl_scan_norm.Rdata")

qtl_peaks_norm = find_peaks(DO_qtl_scan_normal, cross_basic$pmap, threshold=4, drop=1.5)

DO_qtl_scan_binary_norm = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)


qtl_peaks_bin_norm = find_peaks(DO_qtl_scan_binary_norm, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_both_norm = rbind(qtl_peaks_norm,qtl_peaks_bin_norm)

##look at work post yield and toughness. Chr2 , goes away when normalized. Due to extreme phenotypes?
WPY = scan1coef(apr[,2], pheno = cross_basic$pheno[,"bending_work_post_yield"],kinship = k_loco[["2"]], covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
WPY = clean_scan1(WPY)

plot_coefCC(WPY, cross_basic$pmap)

WPY_blup = scan1blup(apr[,2], pheno = cross_basic$pheno[,"bending_work_post_yield"],kinship = k_loco[["2"]], covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)

WPY_blup = clean_scan1(WPY_blup)
plot_coefCC(WPY_blup, cross_basic$pmap)


m = maxmarg(pr,minprob = 0.5,chr = 2, pos=90.856352, map = cross_basic$pmap,return_char = TRUE)
plot_pxg(geno = m, pheno = cross_basic$pheno[,"bending_work_post_yield"],sort = TRUE, force_labels = TRUE)

x = find_marker(cross_basic$pmap,chr = 2,pos = 90.856352)
which(pheno_wpy == sort(pheno_wpy,decreasing = T)[1])
pr$`2`["54",,x]
pr$`2`["341",,x]
#####
#WHAT DO ZERO VALUES MEAN?
pheno_wpy = cross_basic$pheno[,"bending_work_post_yield"]
pheno_wpy = na.omit(pheno_wpy)
pheno_wpy = pheno_wpy[-which(pheno_wpy == sort(pheno_wpy,decreasing = T)[1])]
WPY_q = scan1(apr[,2], pheno = pheno_wpy,kinship = k_loco[["2"]], covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
plot(WPY_q, cross_basic$pmap)
#####
#perms
norm_perm = scan1perm(apr,pheno_combined,k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],n_perm = 1000,perm_Xsp = TRUE,cores = 0)
save(norm_perm,file = "./results/Rdata/norm_perm.Rdata")

perm = scan1perm(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],n_perm = 1000,perm_Xsp = TRUE, cores=0)
save(perm,file = "./results/Rdata/perm.Rdata")


#54 in the middle for slc39a13, near the highest for ptprj, near lowest for Kbtbd4, higher for Celf1, 
#second lowest for ddb2 (lowest is 397, 3rd is 96), 
k = calc_kinship(probs = pr)
herit = est_herit(pheno = pheno_combined,kinship = k,addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])


