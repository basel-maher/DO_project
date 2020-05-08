##supp tables

#S1
#from est_herit


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

load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")

full = combat_annot[,-c(1,3,44)]
full = full[,c(2,1,3:ncol(full))]

colnames(full)[2] = "module"
colnames(full)[2:ncol(full)] = paste0(colnames(full)[2:ncol(full)], "_C")


female = combat_annot_f[,-c(1,3,50)]
female = female[,c(2,1,3:ncol(female))]
colnames(female)[2] = "module"
colnames(female)[2:ncol(female)] = paste0(colnames(female)[2:ncol(female)], "_F")


male = combat_annot_m[,-c(1,3,45)]
male = male[,c(2,1,3:ncol(male))]
colnames(male)[2] = "module"

colnames(male)[2:ncol(male)] = paste0(colnames(male)[2:ncol(male)], "_M")


S4 = merge(full, female,by="Gene.Name", all=T)

S4 = merge(S4, male,by="Gene.Name", all=T)


S4 = S4[,c(1,2,42,88, 3:41, 89:128)]
write.csv(S4, file = "~/Desktop/supp_tables/S4.csv", row.names = F)




#S5 - top 100 GO terms for each network, sorted by pval
load("./results/Rdata/networks/GO_sft4.RData")

networks = as.data.frame(matrix(nrow=100, ncol = 10000))
net_counter = 1
col_counter = 1
for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  networks[,col_counter] = n$Term[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_C_GO_term")
  
  col_counter = col_counter + 1

  networks[,col_counter] = n$classic[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_C_p_value")
  
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
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_M_GO_term")
  
  col_counter = col_counter + 1
  
  networks[,col_counter] = n$classic[1:100]
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_M_p_value")
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}


networks = networks[,c(1:248)]

write.csv(networks, file = "~/Desktop/supp_tables/S5.csv", row.names = F)


#S6
library(Hmisc)
superset = read.delim("./results/flat/superduperset_sansGWAS.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]
superset = capitalize(superset)
write.csv(superset, file = "~/Desktop/supp_tables/S6.csv", row.names = F)




# full_net = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv", stringsAsFactors = FALSE)
# female_net = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv", stringsAsFactors = FALSE)
# male_net = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv", stringsAsFactors = FALSE)
# 
# merge(full_net, female_net, male_net, )


#S7
full = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv", stringsAsFactors = F)
full = full[,c(1:4,16,17)]

colnames(full)[2:ncol(full)] = paste0(colnames(full)[2:ncol(full)], "_C")



female = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv", stringsAsFactors = F)
female = female[,c(1:4,17,18)]

colnames(female)[2:ncol(female)] = paste0(colnames(female)[2:ncol(female)], "_F")


male = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv", stringsAsFactors = F)
male = male[,c(1:4,17,18)]

colnames(male)[2:ncol(male)] = paste0(colnames(male)[2:ncol(male)], "_M")

S7 = merge(full, female, by="gene", all=T)
S7 = merge(S7, male, by="gene", all=T)

S7 = S7[,c(1,2,7,12,3:6,8:11,13:16)]
colnames(S7) = c("gene","module_C","module_F","module_M","num_neib_C","num_bone_neib_C","nominal_pval_C", "FDR_pval_C", "num_neib_F","num_bone_neib_F","nominal_pval_F","FDR_pval_F","num_neib_M","num_bone_neib_M","nominal_pval_M","FDR_pval_M")

write.csv(S7, file = "~/Desktop/supp_tables/S7.csv", row.names = F)



#S8 genes that define rasd1 cluster
library(Seurat)
load("./results/Rdata/seurat_ob.Rdata") #loads as "ob



ob.markers <- FindAllMarkers(ob, only.pos = T)

#Rasd1 cluster (Cluster 10)

ob.markers[which(ob.markers$gene =="Sertad4"),] # 10

#cluster 10 genes
cluster10.markers <- FindMarkers(ob, ident.1 = 10, min.pct = 0.25,only.pos = T)
S8 = cluster10.markers[which(cluster10.markers$p_val_adj <= 0.05),]

write.csv(S8, file = "~/Desktop/supp_tables/S8.csv", row.names = T)


#S9
qtl_loc = read.csv("./results/flat/qtl_loc", stringsAsFactors = FALSE)

#convert to BED3 format while getting max interval size in each CI
bed = as.data.frame(matrix(ncol=3, nrow=length(unique(qtl_loc$locus))))
colnames(bed) = c("chr","start","end")

for(i in 1:length(unique(qtl_loc$locus))){
  sub = subset(qtl_loc, locus == unique(qtl_loc$locus)[i])
  start = min(sub$ci_lo)*1000000
  end = max(sub$ci_hi)*1000000
  chrom = paste0("chr",unique(sub$chr))
  bed[i,] = c(chrom, start, end)
}

lifted = read.table("./results/flat/lifted_qtl_loci.bed")

###
mouse = bed
mouse$region = paste0(mouse$chr,":", mouse$start, "-", mouse$end)

human = lifted
human$region = paste0(human$V1,":", human$V2, "-", human$V3)

S9 = cbind(mouse$region, human$region)

colnames(S9) = c("mouse loci", "syntenic human loci")
write.csv(S9, file = "~/Desktop/supp_tables/S9.csv", row.names = F)


#S10
#SIFT annotations






#S11
#all local eqtl

load("./results/Rdata/local_eqtl.Rdata")
local_eqtl = local_eqtl[,c(1,8,3,4,5,6,7,9,10,11,12,16)]

head(local_eqtl)
local_eqtl$Start = local_eqtl$Start/1000000
local_eqtl$End = local_eqtl$End/1000000

write.csv(local_eqtl, file = "~/Desktop/supp_tables/S11.csv", row.names = F)

#S12
#qsox1 cluster

load("./results/Rdata/seurat_ob.Rdata") #loads as "ob



ob.markers <- FindAllMarkers(ob, only.pos = TRUE)

#Rasd1 cluster (Cluster 10)

ob.markers[which(ob.markers$gene =="Qsox1"),] # 11

# cluster 1 genes
cluster1.markers <- FindMarkers(ob, ident.1 = 1, min.pct = 0.25,only.pos = T)
S12 = cluster1.markers[which(cluster1.markers$p_val_adj <= 0.05),]

write.csv(S12, file = "~/Desktop/supp_tables/S12.csv", row.names = T)


#S13 nad S14 from Larry