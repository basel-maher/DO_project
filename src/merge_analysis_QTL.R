library(dplyr)
##merge analysis
##Do merge analysis for each of the QTL that pass threshold, in the confidence interval for each locus

#DOWNLOAD THE FILES. DOWNLOADED NOV. 14 2019
#download.file("https://ndownloader.figshare.com/files/18533342", "./data/CCdb/cc_variants.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609252", "./data/CCdb/mouse_genes_mgi.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609261", "./data/CCdb/mouse_genes.sqlite")
#

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

#load qtl_loc
qtl_loc = read.csv("./results/flat/qtl_loc")

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

####
#### Normalization####
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

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}

#
#
#


merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")



for(i in 1:nrow(qtl_loc)){
  print(i)
  merge[[i]] = list()
  pheno = qtl_loc$lodcolumn[i]
  chr = qtl_loc$chr[i]
  start = qtl_loc$ci_lo[i]
  end = qtl_loc$ci_hi[i]
  #use same covars as qtl mapping.sex,age,BW and gen
  out_snps <- scan1snps(pr, cross_basic$pmap, pheno_combined[,pheno], k_loco[[chr]],  addcovar =  new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                        query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
  
}

save(merge,file = "./results/Rdata/merge_QTL.Rdata")




#for each merge analysis, take the snps that are within 15% LODs of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)


save(merge_top,file = "./results/Rdata/merge_top_QTL.Rdata")
###
###
###
#Ma.Ar chromosome 3 merge analysis
##
##
##
##
##
#chr3 Ma.Ar QTL Trace
# 
# load("./results/Rdata/DO_qtl_scan_norm.Rdata")
# 
# 
# start = 66.56421 - 0.250
# end = 69.96983 + 0.250
# chr = 3
# out_snps_maar_3 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ma.Ar"], k_loco[["3"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
#                              query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=T)
# 
# 
# variants_locus = query_variants(chr, start, end)
# genes_locus <- query_genes(chr, start, end)
# 
# 
# if("pseudogene" %in% genes_locus$mgi_type){
#   genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
# }
# 
# if("miRNA gene" %in% genes_locus$mgi_type){
#   genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
# }
# 
# if("rRNA gene" %in% genes_locus$mgi_type){
#   genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
# }
# 
# top_maar3 <- top_snps(out_snps_maar_3$lod, out_snps_maar_3$snpinfo, drop = 0.15 * max(out_snps_maar_3$lod))
# top_maar3[order(top_maar3$lod,decreasing = T),]
# 
# 
# #get the colocalizing eqtl
# load("./results/Rdata/local_eqtl.Rdata")#eqtl
# 
# x = local_eqtl[which(local_eqtl$Gene.Name %in% genes_locus$Name),"lodcolumn"]
# x = paste0(x,"_3")
# 
# load("./results/Rdata/merge_top_local_eqtl.Rdata") #merge frame
# 
# merge_top_maar = merge_top[which(names(merge_top) %in% x)]
# 
# coloc_snps = c()
# for(i in 1:length(merge_top_maar)){
#   if(any(top_maar3$snp_id %in% merge_top_maar[[i]]$snp_id)){
#     print(names(merge_top_maar)[i])
#     coloc_snps = append(coloc_snps,top_maar3[which(top_maar3$snp_id %in% merge_top_maar[[i]]$snp_id),"snp_id"])
#   }
# }
# 
# 
# coloc_snps = unique(coloc_snps) #colocalizing snps
# 
# 
