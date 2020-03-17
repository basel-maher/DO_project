library(ggplot2)
library(ggrepel)
library(ggthemes)
library(qtl2)
load("./data/POMP/model.probs1_qtl2_cleaned.Rdata")
load("./data/POMP/Xcovar.Rdata")

load("./data/POMP/MM_snps1_pmap.Rdata")


pheno_pomp = read.csv("./data/POMP/pheno1.csv", stringsAsFactors = FALSE)

apr_pomp = genoprob_to_alleleprob(model.probs1_qtl2_cleaned)

k_loco_pomp <- calc_kinship(apr_pomp, "loco")
rownames(pheno_pomp) = pheno_pomp$X

scan1_pomp = scan1(apr_pomp, pheno_pomp[,c(11:17)], k_loco_pomp, Xcovar=Xcovar, addcovar = pheno_pomp[,c("Sex", "Diet","AgeSac","MassSac")])

pomp_peaks = find_peaks(scan1_pomp, MM_snps1_pmap, threshold=4, drop=1.5)

plot(scan1_pomp, MM_snps1_pmap$`1`,lodcolumn = 2)

abline(h = 7.8, col = "red", lwd = 2)

ML_pomp = pheno_pomp$ML
names(ML_pomp) = pheno_pomp$X

coef_ML_pomp <- scan1coef(apr_pomp[,"1"], ML_pomp, kinship = k_loco_pomp[["1"]], addcovar = pheno_pomp[,c("Sex", "Diet","AgeSac","MassSac")])
plot_coefCC(coef_ML_pomp, MM_snps1_pmap["1"], scan1_output=subset(scan1_pomp, lodcolumn=2))

coef_ML_pomp_blup = scan1blup(apr_pomp[,"1"], ML_pomp, kinship = k_loco_pomp[["1"]], addcovar = pheno_pomp[,c("Sex", "Diet","AgeSac","MassSac")])
plot_coefCC(coef_ML_pomp_blup, MM_snps1_pmap["1"], scan1_output=subset(scan1_pomp, lodcolumn=2),legend = "topleft")

