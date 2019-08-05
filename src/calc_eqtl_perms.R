##calculate permutations for eQTL analysis.
#because the counts are qnormed, they have the same distribution. 
#This allows us to sample a small population of them, permute them, and use LOD threshold from those as our global lod thresholds

library(qtl2)
#for mapping, you need the allele probs, kinship and cross file
load("./results/Rdata/apr_basic_cleaned.Rdata")
load("./results/Rdata/paper_analyses/cross_eqtl.Rdata")
load("./results/Rdata/k_loco_basic_cleaned.Rdata")

#randomly choose 50 genes to permute and save them
#perm_cols = sample(1:ncol(cross_eqtl$pheno),size = 50, replace = FALSE)
#saveRDS(perm_cols, "./results/Rdata/perm_cols.RDS")
cols = readRDS("./results/Rdata/perm_cols.RDS")


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
#eQTL mapping - local perms

set.seed(8675309)

#use sex and first 35 PEER covars
for(num in cols){
  name = colnames(cross_eqtl$pheno)[num]
  
  perm = scan1perm(apr, cross_eqtl$pheno[,num],k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,17:51)],n_perm = 1000,perm_Xsp = TRUE,chr_lengths=chr_lengths(cross_eqtl$gmap),cores=20)

  save(perm,file = as.character(paste0("./results/Rdata/eqtl_perms/localeqtl_perm_",name,".Rdata")))
}
##

#eQTL mapping - distal perms

set.seed(8675309)

#use all 48 peer covars
for(num in cols){
  name = colnames(cross_eqtl$pheno)[num]
  
  perm = scan1perm(apr, cross_eqtl$pheno[,num],k_loco, Xcovar=Xcovar, addcovar = covar[,c(17:ncol(covar))],n_perm = 1000,perm_Xsp = TRUE,chr_lengths=chr_lengths(cross_eqtl$gmap),cores=20)
  
  save(perm,file = as.character(paste0("./results/Rdata/eqtl_perms/distaleqtl_perm_",name,".Rdata")))
}
##

#read the output and get a permutation value. put them into a dataframe
perm_frame = as.data.frame(matrix(ncol = 5, nrow=50))

for(i in 1:length(cols)){
  name = colnames(cross_eqtl$pheno)[cols[i]]
  load(paste0("./results/Rdata/eqtl_perms/localeqtl_perm_",name,".Rdata"))#loads as "perm"
  autosome_lod = summary(perm, 0.05)$A[[1]] #autosomal lod threshold, alpha=0.05
  x_lod = summary(perm, 0.05)$X[[1]] #X chrom lod threshold, alpha=0.05
  perm_frame[i,1] = name
  perm_frame[i,2] = autosome_lod
  perm_frame[i,3] = x_lod
  
  load(paste0("./results/Rdata/eqtl_perms/distaleqtl_perm_",name,".Rdata"))#do the same for distal eqtl
  autosome_lod = summary(perm, 0.05)$A[[1]]
  x_lod = summary(perm, 0.05)$X[[1]]
  
  perm_frame[i,4] = autosome_lod
  perm_frame[i,5] = x_lod
}
colnames(perm_frame) = c("gene","local_autosomal", "local_X","distal_autosomal","distal_X")
