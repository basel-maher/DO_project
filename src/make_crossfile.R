##create cross file for qtl2 analysis
library(qtl2)
load("./results/Rdata/full_pheno_table.Rdata")
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata")) #GigaMUGA SNPs
#broman
dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/master/UWisc/"
gm <- read.csv(paste0(dir, "gm_uwisc_v1.csv"))
#remove the SNPs that dont map uniquely to the genome (broman)
GM_snps = GM_snps[-match(gm[which(gm$unique==FALSE),"marker"],GM_snps$marker),]
####

geno = readRDS("./results/GIGAMUGA/geno.final_merged.RDS")

#GENERATE APPROPRIATE PHENO FILES
do_pheno = full_pheno_table
#remove categorical vars
do_pheno = do_pheno[,-which(names(do_pheno) %in% c("sac_date","coat_color","DOB","comments","DO_generation"))]

#which colnames have periods in the geno files? (These have periods because they were re-genotyped)
colnames(geno)[grep(x = colnames(geno),pattern = "\\.")]
#add periods to the pheno file so it matches the geno name
do_pheno$Mouse.ID[c(265,271,291,324,346,350,352,4,371)] = c("265.1","271.1","291.1","324.1","346.1","350.1","352.1","4.1","371.1")


#make sure all numeric and rownames are correct
rname_pheno = rownames(do_pheno) = do_pheno$Mouse.ID
do_pheno = as.data.frame(lapply(do_pheno,as.numeric))
rownames(do_pheno) = do_pheno$Mouse.ID = rname_pheno
 
#subset only those that have genotyping data
do_pheno = do_pheno[which(do_pheno$Mouse.ID %in% colnames(geno)),]
write.csv(do_pheno, "./results/flat/do_pheno.csv", quote = F)
###
###

#prepare covariate data. Must contain sex column and a column giving the DO generation
ix = which(colnames(full_pheno_table) %in% c("Mouse.ID","sex","sac_date","age_at_sac_days","body_weight","body_length","DO_generation"))

#actually generation needs to be a "surrogate variable". 0s and 1s, one less variable than the number of generations
generation <- factor(full_pheno_table$DO_generation)
X <- model.matrix(~ generation)[,-1]

do_covar = full_pheno_table[,ix]
do_covar = cbind(do_covar,X)#add generation surrogate vars

#same subtitution of rownames as above
colnames(geno)[grep(x = colnames(geno),pattern = "\\.")]
do_covar$Mouse.ID[c(265,271,291,324,346,350,352,4,371)] = c("265.1","271.1","291.1","324.1","346.1","350.1","352.1","4.1","371.1")

colnames(do_covar)[7] = "ngen"#rename
do_covar$ngen = apply(do_covar, 1, function(x) strsplit(x[7], split = "G")[[1]][2])#remove "G

rownames(do_covar) = do_covar$Mouse.ID 

#use only those with genotyping info
do_covar = do_covar[which(do_covar$Mouse.ID %in% colnames(geno)),]
write.csv(do_covar, "./results/flat/do_covar.csv", quote = FALSE,row.names = FALSE)



#prepare control file. control file paths relative to control file path
setwd("./results/flat/")

chr <- c(1:19, "X")
write_control_file("./control_file_basic.json",
                   crosstype="do",
                   description="control file for QTL analysis, basic",
                   founder_geno_file=paste0("../../data/GIGAMUGA/GM_processed_files/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("../../data/GIGAMUGA/GM_processed_files/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("../../data/GIGAMUGA/GM_processed_files/GM_pmap", chr, ".csv"),
                   geno_file=paste0("../GIGAMUGA/qtl2_batches_1-4/DO_qtl2_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="../flat/do_pheno.csv",
                   covar_file="../flat/do_covar.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen",
                   na.strings=list(NA,"-"),overwrite = T)
#


#read the cross file
#112869 markers
cross = read_cross2("./control_file_basic.json")



#get the appropriate marker set for genotype probability calculation
GM_snps_tier1_2 = subset(GM_snps, tier %in% c(1,2))


cross_combined = pull_markers(cross, GM_snps_tier1_2$marker)

#109910 markers remaining - final
cross_combined = pull_markers(cross_combined, rownames(geno.final_merged))

cross_basic = cross_combined
save(cross_basic, file ="../../results/Rdata/cross_basic.Rdata")

##

