library(ggplot2)
library(ggrepel)
library(ggthemes)
library(qtl2)
load("./data/POMP/model.probs1_qtl2_cleaned.Rdata")
load("./data/POMP/Xcovar.Rdata")

load("./data/POMP/MM_snps1_pmap.Rdata")

do_pmap = cross_basic$pmap
##
for(lst in 1:length(MM_snps1_pmap)){
  markers = names(MM_snps1_pmap[[lst]][which(names(MM_snps1_pmap[[lst]]) %in% names(do_pmap[[lst]]))])
  MM_snps1_pmap[[lst]][markers] == do_pmap[[lst]][markers]
  
  if(all(MM_snps1_pmap[[lst]][markers] == do_pmap[[lst]][markers])){
    print(lst)
    do_pmap[[lst]] == MM_snps1_pmap[[lst]]
  }
}
##one is different, on chromosome 10 (UNC17323847)
#everything else is the same btwn DO and POMP maps. so set all to be the same and remove UNC17323847 

##

pheno_pomp = read.csv("./data/POMP/pheno1.csv", stringsAsFactors = FALSE)

apr_pomp = genoprob_to_alleleprob(model.probs1_qtl2_cleaned)

k_loco_pomp <- calc_kinship(apr_pomp, "loco")
rownames(pheno_pomp) = pheno_pomp$X

scan1_pomp = scan1(apr_pomp, pheno_pomp[,c(11:17)], k_loco_pomp, Xcovar=Xcovar, addcovar = pheno_pomp[,c("Sex", "Diet","AgeSac","MassSac")])
save(scan1_pomp, file="./results/Rdata/scan1_pomp.Rdata")
pomp_peaks = find_peaks(scan1_pomp, do_pmap, threshold=4, drop=1.5)

plot(scan1_pomp, MM_snps1_pmap$`1`,lodcolumn = 2)

abline(h = 7.8, col = "red", lwd = 2)

ML_pomp = pheno_pomp$ML
names(ML_pomp) = pheno_pomp$X

coef_ML_pomp <- scan1coef(apr_pomp[,"1"], ML_pomp, kinship = k_loco_pomp[["1"]], addcovar = pheno_pomp[,c("Sex", "Diet","AgeSac","MassSac")])
plot_coefCC(coef_ML_pomp, do_pmap["1"], scan1_output=subset(scan1_pomp, lodcolumn=2))

coef_ML_pomp_blup = scan1blup(apr_pomp[,"1"], ML_pomp, kinship = k_loco_pomp[["1"]], addcovar = pheno_pomp[,c("Sex", "Diet","AgeSac","MassSac")])
plot_coefCC(coef_ML_pomp_blup, do_pmap["1"], scan1_output=subset(scan1_pomp, lodcolumn=2),legend = "topleft")
save(coef_ML_pomp_blup, file="./results/Rdata/coef_ML_pomp_blup.Rdata")

