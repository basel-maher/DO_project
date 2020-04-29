##supp tables

#S1



#S2
#Holm-Bonferroni correlations with femoral strength
library(psych)
set.seed(8675309)
#lasso regression to find best predictors of max load
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")


#normalize
norm_pheno = as.data.frame(cross_basic$pheno)

norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1 + 1
norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2 + 1
norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1

norm_pheno$bending_work_post_yield = norm_pheno$bending_work_post_yield + 1
norm_pheno$bending_disp_at_yield = norm_pheno$bending_disp_at_yield + 1
norm_pheno$bending_PYD = norm_pheno$bending_PYD + 1

norm_pheno = as.data.frame(log10(norm_pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,52,54:58,61:70,72,74,76)]))

pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,36,42,50,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}

##########


#remove some traits
pheno_combined = pheno_combined[,-c(1:9,54:57)]

pheno_combined=pheno_combined[,c(1,2,45,3:10,46:49, 11:26, 50,51,41:44,27:40,52:55)]
#pearson cor . adjusted for multiple comparisons (Holm)
x = psych::corr.test(pheno_combined)

#p vals
p = x$p[,"bending_max_load"]
#correlation
r = x$r[,"bending_max_load"]

#S.E.
se = x$se[,"bending_max_load"]
table = cbind(r,p, se)


write.csv(table, file = "~/Desktop/supp_tables/S2.csv")



#S3
table = x$r
p = x$p

for( i in 1:length(table)){
  table[i] = paste0(signif(as.numeric(table[i],3)), "(", signif(as.numeric(p[i],3)), ")")
}

write.csv(table, file = "~/Desktop/supp_tables/S3.csv")


#S4
library(data.table)
library(rowr)
load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")

full = combat_annot[,c("Gene.Name","color")]

full$color = paste0(full$color, "_C")

S4 = data.table::dcast(setDT(full), rowid(color)~color, value.var="Gene.Name")
S4 = S4[,-1]

female = combat_annot_f[,c("Gene.Name","color")]
female$color = paste0(female$color, "_F")

x = data.table::dcast(setDT(female), rowid(color)~color, value.var="Gene.Name")
x = x[,-1]

S4 = cbind.fill(S4, x, MoreArgs = list(fill=NA))
S4 = S4[,-85] #"object"


male = combat_annot_m[,c("Gene.Name","color")]
male$color = paste0(male$color, "_M")

x = data.table::dcast(setDT(male), rowid(color)~color, value.var="Gene.Name")
x = x[,-1]

S4 = cbind.fill(S4, x, MoreArgs = list(fill=NA))
S4 = S4[,-125] #"object"

write.csv(S4, file = "~/Desktop/supp_tables/S4.csv")




#S5 - top 100 GO terms for each network, sorted by pval
load("./results/Rdata/networks/GO_sft4.RData")

networks = as.data.frame(matrix(nrow=100, ncol = 10000))
net_counter = 1
col_counter = 1
for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  networks[,col_counter] = n$Term[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_GO_term")
  
  col_counter = col_counter + 1

  networks[,col_counter] = n$classic[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_p_value")
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}



load("./results/Rdata/networks/GO_Females_sft4.RData")
net_counter = 1

for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  networks[,col_counter] = n$Term[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_F_GO_term")
  
  col_counter = col_counter + 1
  
  networks[,col_counter] = n$classic[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_F_p_value")
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}





load("./results/Rdata/networks/GO_Males_sft5.RData")


net_counter = 1

for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  networks[,col_counter] = n$Term[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_F_GO_term")
  
  col_counter = col_counter + 1
  
  networks[,col_counter] = n$classic[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_F_p_value")
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}


networks = networks[,c(1:248)]

write.csv(networks, file = "~/Desktop/supp_tables/S5.csv", row.names = F)






# full_net = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv", stringsAsFactors = FALSE)
# female_net = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv", stringsAsFactors = FALSE)
# male_net = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv", stringsAsFactors = FALSE)
# 
# merge(full_net, female_net, male_net, )
