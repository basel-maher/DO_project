##create cross file for qtl2 analysis
library(qtl2)
load("./results/Rdata/full_pheno_table.Rdata")
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
 

write.csv(do_pheno, "./results/flat/do_pheno.csv", quote = F)
###
###

#prepare covariate data. Must contain sex dolumn and a column giving the DO generation
ix = which(colnames(full_pheno_table_ordered) %in% c("Mouse.ID","sex","Generation","Sac.date","Age.at.Sac..days.","Body.Weight","Nose.Anal.Length"))

do_covar = full_pheno_table_ordered[,ix]

colnames(geno.final_merged)[grep(x = colnames(geno.final_merged),pattern = "\\.")]
do_covar$Mouse.ID[c(265,271,291,324,346,350,352,4)] = c("265.1","271.1","291.1","324.1","346.1","350.1","352.1","4.1")

colnames(do_covar)[7] = "ngen"
do_covar$ngen = apply(do_covar, 1, function(x) strsplit(x[7], split = "G")[[1]][2])

do_covar = do_covar[sort(c(371,which(do_covar$Mouse.ID %in% colnames(geno.final_merged)))),]#add 371 because its good now (2019 batch fixed dec. 2018)

rownames(do_covar) = do_covar$Mouse.ID 
write.csv(do_covar, "~/Desktop/DO_proj/data/GIGAMUGA/do_covar.csv", quote = FALSE,row.names = FALSE)



#prepare control file. control file paths relative to control file path
chr <- c(1:19, "X")
write_control_file("~/Desktop/DO_proj/data/GIGAMUGA/control_file_gigamuga.json",
                   crosstype="do",
                   description="control file for eQTL analysis",
                   founder_geno_file=paste0("GM_processed_files/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM_processed_files/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM_processed_files/GM_pmap", chr, ".csv"),
                   geno_file=paste0("qtl2_batches_1-2-3-4/qtl2_batches_1-2-3-4_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="do_pheno.csv",
                   covar_file="do_covar.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen",
                   na.strings=list(NA,"-"),overwrite = T)
#
#read the cross file
#114184 markers
cross = read_cross2("~/Desktop/DO_proj/data/GIGAMUGA/control_file_gigamuga.json")



#get the appropriate marker set for genotype probability calculation
GM_snps_tier1_2 = subset(GM_snps, tier %in% c(1,2))


#112400 markers
cross_combined = pull_markers(cross, GM_snps_tier1_2$marker)
#geno.final_merged from array_QC.R
#110545 markers remaining - final
cross_combined = pull_markers(cross_combined, rownames(geno.final_merged))

save(cross_combined, file ="~/Desktop/DO_proj/data/GIGAMUGA/cross_combined.Rdata")

##

#load("~/Desktop/DO_proj/data/GIGAMUGA/cross_combined.Rdata")
#####
####

covar = as.matrix(cross_combined$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1
covar = covar[,c(2,3,4,5,6)]
covar = apply(covar,2,as.numeric)
#covar[,c(3,4,5)] = apply(covar[,c(3,4,5)],2,as.numeric)
rownames(covar) = rownames(cross_combined$covar)
#####

