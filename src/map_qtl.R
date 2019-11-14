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
}
#only 15,18,19,20,22,33,36,42,50,53,59 are Normal



#qtl mapping for continuous traits

#which traits require FL as covar? bending traits?
DO_qtl_scan = scan1(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
#save(DO_qtl_scan,file = "./results/Rdata/DO_qtl_scan.Rdata")


DO_qtl_scan_binary = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
#save(DO_qtl_scan_binary,file = "./results/Rdata/DO_qtl_scan_binary.Rdata")
load("./results/Rdata/DO_qtl_scan_binary.Rdata")

#find peaks and bind them together
qtl_peaks = find_peaks(DO_qtl_scan, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_binary = find_peaks(DO_qtl_scan_binary, cross_basic$pmap, threshold=4, drop=1.5)

qtl_peaks_both = rbind(qtl_peaks,qtl_peaks_binary)
write.csv(qtl_peaks_both, file = "./results/flat/qtl_peaks.csv",row.names = FALSE,quote = FALSE)


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

####remove glucose measurements for decap, non fasted, etc

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



#########
##look at work post yield and toughness. Chr2 , goes away when normalized. Due to extreme phenotypes?
# WPY = scan1coef(apr[,2], pheno = cross_basic$pheno[,"bending_work_post_yield"],kinship = k_loco[["2"]], covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# WPY = clean_scan1(WPY)
# 
# plot_coefCC(WPY, cross_basic$pmap)
# 
# WPY_blup = scan1blup(apr[,2], pheno = cross_basic$pheno[,"bending_work_post_yield"],kinship = k_loco[["2"]], covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# 
# WPY_blup = clean_scan1(WPY_blup)
# plot_coefCC(WPY_blup, cross_basic$pmap)
# 
# 
# m = maxmarg(pr,minprob = 0.5,chr = 2, pos=90.856352, map = cross_basic$pmap,return_char = TRUE)
# plot_pxg(geno = m, pheno = cross_basic$pheno[,"bending_work_post_yield"],sort = TRUE, force_labels = TRUE)
# 
# x = find_marker(cross_basic$pmap,chr = 2,pos = 90.856352)
# which(pheno_wpy == sort(pheno_wpy,decreasing = T)[1])
# pr$`2`["54",,x]
# pr$`2`["341",,x]
# #####
# #WHAT DO ZERO VALUES MEAN?
# pheno_wpy = cross_basic$pheno[,"bending_work_post_yield"]
# pheno_wpy = na.omit(pheno_wpy)
# pheno_wpy = pheno_wpy[-which(pheno_wpy == sort(pheno_wpy,decreasing = T)[1])]
# WPY_q = scan1(apr[,2], pheno = pheno_wpy,kinship = k_loco[["2"]], covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
# plot(WPY_q, cross_basic$pmap)
# #####
# #perms
# norm_perm = scan1perm(apr,pheno_combined,k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],n_perm = 1000,perm_Xsp = TRUE,cores = 0)
# save(norm_perm,file = "./results/Rdata/norm_perm.Rdata")
# 
# perm = scan1perm(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],n_perm = 1000,perm_Xsp = TRUE, cores=0)
# save(perm,file = "./results/Rdata/perm.Rdata")
# 
# 
# #54 in the middle for slc39a13, near the highest for ptprj, near lowest for Kbtbd4, higher for Celf1, 
# #second lowest for ddb2 (lowest is 397, 3rd is 96), 
# k = calc_kinship(probs = pr)
# herit = est_herit(pheno = pheno_combined,kinship = k,addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])
# 
# ###
# #remove NAs from pheno, then remove empty covars
# 
####
### interactive model w/ sex
DO_qtl_int = scan1(apr, cross_basic$pheno[,c(6:70,72,74,76)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = covar[,"sex"],cores = 2)
#save(DO_qtl_int,file = "./results/Rdata/DO_qtl_int.Rdata")
load("./results/Rdata/DO_qtl_int.Rdata")
qtl_peaks_int = find_peaks(DO_qtl_int, cross_basic$pmap, threshold=4, drop=1.5)


DO_qtl_int_binary = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = covar[,"sex"],cores = 2)
#save(DO_qtl_int_binary,file = "./results/Rdata/DO_qtl_int_binary.Rdata")
load("./results/Rdata/DO_qtl_int_binary.Rdata")

qtl_peaks_binary_int = find_peaks(DO_qtl_int_binary, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_both_int = rbind(qtl_peaks_int,qtl_peaks_binary_int)


##diff both
x = merge(qtl_peaks_both_int,qtl_peaks_both_norm, by=c("lodcolumn","chr"))

#same but norm
DO_qtl_int_norm = scan1(apr, pheno_combined, k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = new_covar[,"sex"],cores = 2)
save(DO_qtl_int_norm,file = "./results/Rdata/DO_qtl_int_norm.Rdata")
load("./results/Rdata/DO_qtl_int_norm.Rdata")

qtl_peaks_int_norm = find_peaks(DO_qtl_int_norm, cross_basic$pmap, threshold=4, drop=1.5)

DO_qtl_int_binary_norm = scan1(apr, cross_basic$pheno[,c(71,73,75,77)], Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],intcovar = new_covar[,"sex"],cores = 2)
save(DO_qtl_int_binary_norm,file = "./results/Rdata/DO_qtl_int_binary_norm.Rdata")
load("./results/Rdata/DO_qtl_int_binary_norm.Rdata")
qtl_peaks_int_bin_norm = find_peaks(DO_qtl_int_binary_norm, cross_basic$pmap, threshold=4, drop=1.5)

qtl_peaks_both_int = rbind(qtl_peaks_int_bin_norm,qtl_peaks_int_norm)

#Do LODi = LODf-LODa. F is full model (additive + int), a is additive only





#####
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

##same for norm
qtl_peaks_both_norm$perm_thresh = NA
perm = list.files("./results/Rdata/qtl_perms/")

perm_files = perm[grep("norm_perms_",perm)]

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
qtl = qtl_peaks_both[which(qtl_peaks_both$lod >= qtl_peaks_both$perm_thresh),]
qtl_norm = qtl_peaks_both_norm[which(qtl_peaks_both_norm$lod >= qtl_peaks_both_norm$perm_thresh),]
###

TMD_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Ct.TMD"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(TMD_blup, file = "./results/Rdata/TMD_blup.Rdata")

#load("./results/Rdata/TMD_blup.Rdata")
png("~/Desktop/TMD.png")
plot_coefCC(TMD_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=45), main = "TMD")
dev.off() 


plot_coefCC(TMD_blup, cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=45), main = "TMD")

TMD_coef_scan = scan1coef(apr,pheno = pheno_combined[,"uCT_Ct.TMD"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)

pMOI_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_pMOI"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(pMOI_blup, file = "./results/Rdata/pMOI_blup.Rdata")
load("./results/Rdata/pMOI_blup.Rdata")

png("~/Desktop/pMOI.png")
plot_coefCC(pMOI_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=47), main = "pMOI")
dev.off()

imax_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Imax"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(imax_blup, file = "./results/Rdata/imax_blup.Rdata")
load("./results/Rdata/imax_blup.Rdata")

png("~/Desktop/imax.png")
plot_coefCC(imax_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=48), main = "Imax")
dev.off()

ttAr_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Tt.Ar"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(ttAr_blup, file = "./results/Rdata/ttAr_blup.Rdata")
load("./results/Rdata/ttAr_blup.Rdata")

png("~/Desktop/ttar.png")
plot_coefCC(ttAr_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=43), main = "Tt.Ar")
dev.off()

maAr_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Ma.Ar"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(maAr_blup, file = "./results/Rdata/maAr_blup.Rdata")
load("./results/Rdata/maAr_blup.Rdata")

png("~/Desktop/Ma.Ar.png")
plot_coefCC(maAr_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=42), main = "Ma.Ar")
dev.off()

ML_blup = scan1blup(apr,pheno = pheno_combined[,"ML" ], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(ML_blup, file = "./results/Rdata/ML_blup.Rdata")
load("./results/Rdata/ML_blup.Rdata")

png("~/Desktop/ml.png")
plot_coefCC(ML_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=10),main = "ML")
dev.off()

porosity_blup = scan1blup(apr,pheno = pheno_combined[,"uCT_Ct.porosity"], kinship = k_loco[["1"]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
save(porosity_blup, file = "./results/Rdata/porosity_blup.Rdata")
load("./results/Rdata/porosity_blup.Rdata")

png("~/Desktop/porosity.png")
plot_coefCC(porosity_blup[4000:7500,], cross_basic$pmap,legend = "topleft",scan1_output = subset(DO_qtl_scan_normal, lodcolumn=46),main="Ct.Porosity")
dev.off()

#
TMD_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
TMD_coef = TMD_blup[TMD_mark,1:8]

pMOI_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
pMOI_coef = pMOI_blup[pMOI_mark,1:8]

imax_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
imax_coef = imax_blup[imax_mark,1:8]

ttAr_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.19293)
ttAr_coef = ttAr_blup[ttAr_mark,1:8]

maAr_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.32971)
maAr_coef = maAr_blup[maAr_mark,1:8]


ML_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.35757)
ML_coef = ML_blup[ML_mark,1:8]

por_mark = find_marker(cross_basic$pmap, chr = 1, pos = 156.18734)
por_coef = porosity_blup[por_mark,1:8]

cor(pMOI_coef, TMD_coef, method = "k")

###


load("./results/Rdata/cross_eqtl.Rdata")



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

qsox1_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1313"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
save(qsox1_blup, file = "./results/Rdata/qsox1_blup.Rdata")

png("~/Desktop/qsox1.png")
plot_coefCC(qsox1_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Qsox1")
dev.off()
#load("./results/Rdata/TMD_blup.Rdata")
#Lhx4 :NA
#Acbd6: MSTRG.1312
Acbd6_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1312"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/acbd6.png")
plot_coefCC(Acbd6_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Acbd6")
dev.off()
#Cep350: MSTRG.1316	

Cep350_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1316"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/cep350.png")
plot_coefCC(Cep350_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Cep350")
dev.off()
#$Gm37571: NA
#Gm37539: NA 
#Xpr1: MSTRG.1309	
Xpr1_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1309"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/xpr1.png")
plot_coefCC(Xpr1_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Xpr1")
dev.off()

#Tor1aip2: MSTRG.1324
Tor1aip2_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1324"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/tor1aip2.png")
plot_coefCC(Tor1aip2_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Tor1aip2")
dev.off()
#Tdrd5: ENSMUSG00000060985
Tdrd5_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"ENSMUSG00000060985"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/tdrd5.png")
plot_coefCC(Tdrd5_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Tdrd5")
dev.off()
#Stx6: MSTRG.1305	
Stx6_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1305"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/stx6.png")
plot_coefCC(Stx6_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Stx6")
dev.off()
#Mr1: MSTRG.1304	
Mr1_blup = scan1blup(apr,pheno = cross_eqtl$pheno[,"MSTRG.1304"], kinship = k_loco[["1"]], addcovar = covar[,c(2,17:51)],cores = 2)
png("~/Desktop/mr1.png")
plot_coefCC(Mr1_blup[4000:7500,], cross_eqtl$pmap,legend = "topleft", main = "Mr1")
dev.off()

##
TMD_mark = find_marker(cross_basic$pmap, chr = 1, pos = 155.10491)
TMD_coef = TMD_blup[TMD_mark,1:8]

mr1_mark = find_marker(cross_eqtl$pmap, chr = 1, pos = 155.10491)
mr1_coef = Mr1_blup[mr1_mark,1:8]

cor(TMD_coef, mr1_coef, method = "k")

