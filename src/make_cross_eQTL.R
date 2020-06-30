#make the cross file for eQTL mapping
library(qtl2)
load("./results/Rdata/cross_basic_cleaned.Rdata") #cross_basic
load("./results/Rdata/full_pheno_table.Rdata") #pheno
load("./results/Rdata/counts_vst_qnorm.Rdata") #normalize counts
geno = readRDS("./results/GIGAMUGA/geno.final_merged.RDS")

#Make the cross file#

#prepare covariate data. Must contain sex column and a column giving the DO generation
#get from previously made qtl cross file
do_covar_eqtl = cross_basic$covar 
#match rows to only get indiviuals with RNAseq data
do_covar_eqtl = do_covar_eqtl[match(rownames(dds.vst),rownames(do_covar_eqtl)),]

#add "Mouse.ID column"
do_covar_eqtl$Mouse.ID = rownames(do_covar_eqtl)
do_covar_eqtl = do_covar_eqtl[,c(17,1:16)]#put ID column first

#read in the PEER factors. These were calculated for 48 factors using no covariates and no intercept 
peer_factors = read.delim("./results/flat/factors_vst_qnorm_redo.txt",stringsAsFactors = FALSE, row.names = rownames(dds.vst), sep=" ", header = FALSE)
for(i in 1:48){colnames(peer_factors)[i] = paste0("PEER_factor_",i)}#rename peer factors
#match the PEER factors to the covar object rownames
peer_factors = peer_factors[match(do_covar_eqtl$Mouse.ID, rownames(peer_factors)),]

#add PEER factors to covar object 
do_covar_eqtl = cbind(do_covar_eqtl,peer_factors)
#remove gens 25,26,30-33, not in data
do_covar_eqtl = do_covar_eqtl[,-c(9,10,14:17)]

#save
write.csv(do_covar_eqtl, "./results/flat/do_covar_eqtl_REDO.csv", quote = FALSE,row.names = FALSE)


#prepare control file. control file paths relative to control file path. In this case, pheno are the normalized counts
setwd("./results/flat/")

chr <- c(1:19, "X")
write_control_file("./control_file_eqtl_REDO.json",
                   crosstype="do",
                   description="control file for eQTL analysis",
                   founder_geno_file=paste0("../../data/GIGAMUGA/GM_processed_files/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("../../data/GIGAMUGA/GM_processed_files/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("../../data/GIGAMUGA/GM_processed_files/GM_pmap", chr, ".csv"),
                   geno_file=paste0("../GIGAMUGA/qtl2_batches_1-4/DO_qtl2_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="../flat/counts_vst_qnorm.csv",
                   covar_file="../flat/do_covar_eqtl_REDO.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen",
                   na.strings=list(NA,"-"),overwrite = T)
#


#read the cross file
cross_eqtl = read_cross2("./control_file_eqtl_REDO.json")

#get same markers as in cross_basic (the QTL cross which was QC'd)
markers = marker_names(cross_basic)

cross_eqtl = pull_markers(cross_eqtl,markers)

save(cross_eqtl, file = "../Rdata/cross_eqtl_REDO.Rdata")

