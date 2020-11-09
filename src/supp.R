library(rtracklayer)
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


S4 = S4[,c(1,2,42,88, 3:41, 43:87,89:128)]

gtf = rtracklayer::import("results/flat/RNA-seq/mus_stringtie_merged.gtf")
gtf = as.data.frame(gtf)
gtf = gtf[,c("gene_name","ref_gene_id")]
gtf = as.data.frame(unique(gtf))
gtf = gtf[-which(is.na(gtf$ref_gene_id)),]
require(dplyr)
gtf.merged = gtf %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ref_gene_id = paste(ref_gene_id, collapse = ","))


S4$name2 = sapply(strsplit(S4$Gene.Name, split = "_isoform"),"[",1)

S4 = merge(S4, gtf.merged, by.x = "name2", by.y = "gene_name")
S4 = S4[-1]
S4 = S4[,c(1,129, 2:128)]
colnames(S4)[2] = "Ensembl_ID"

mart = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations

out <- (getBM(attributes=c('mgi_id','external_gene_name',"external_synonym"),
              filters = 'external_gene_name', values = S4$Gene.Name, mart = mart))


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
#these are mouse genes
library(Hmisc)
superset = read.delim("./results/flat/superduperset_sansGWAS.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]
#superset = capitalize(superset)

mart = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations

out <- (getBM(attributes=c('mgi_id','external_gene_name'),
                            filters = 'external_gene_name', values = superset, mart = mart))

out$external_gene_name = tolower(out$external_gene_name)
superset = as.data.frame(superset)
superset = merge(superset, out, by.x = "superset", by.y="external_gene_name", all.x=T)

mgi = read.delim("./data/MGIhdpQuery_markers_20190728_224719.txt",stringsAsFactors = FALSE)
mgi_mouse = mgi[which(mgi$Organism=="mouse"),]
mgi_mouse$Gene.Symbol = tolower(mgi_mouse$Gene.Symbol)

superset= merge(superset, mgi_mouse, by.x="superset", by.y="Gene.Symbol", all.x=T,)

superset[which(is.na(superset$mgi_id)),"mgi_id"] = superset[which(is.na(superset$mgi_id)),"ID"]

superset[1052,"mgi_id"] = 'MGI:1195971' #pira1
superset[1054,"mgi_id"] = 'MGI:1195974'  #pira6

which(duplicated(superset$superset)) #ctla4 duplicated but no other entries for it
superset = superset[-335,]

superset = superset[,c(1,2)]
colnames(superset) = c("Feature","MGI_ID")
superset$Feature = capitalize(superset$Feature)
##CHECK RIKEN GENE NOMENCLATURE CAPITALIZATION
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

gtf = rtracklayer::import("results/flat/RNA-seq/mus_stringtie_merged.gtf")
gtf = as.data.frame(gtf)
gtf = gtf[,c("gene_name","ref_gene_id")]
gtf = as.data.frame(unique(gtf))
gtf = gtf[-which(is.na(gtf$ref_gene_id)),]
require(dplyr)
gtf.merged = gtf %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ref_gene_id = paste(ref_gene_id, collapse = ","))


S7$name2 = sapply(strsplit(S7$gene, split = "_isoform"),"[",1)

S7 = merge(S7, gtf.merged, by.x = "name2", by.y = "gene_name", all.x=T)

S7 = S7[-1]
S7 = S7[,c(1,17, 2:16)]
colnames(S7)[2] = "Ensembl_ID"

which(is.na(S7$Ensembl_ID))

write.csv(S7, file = "~/Desktop/supp_tables/S7.csv", row.names = F)



#S8 correlation between modules and bone phenotypes
rows = match(colnames(pheno_combined), colnames(moduleTraitCor))

cor = t(moduleTraitCor)[rows,]
pval = t(moduleTraitPvalue)[rows,]

c = cor
for(i in 1:length(c)){
  c[i] = paste0(cor[i], " (", pval[i], ")")
}  

colnames(c) = sapply(strsplit(colnames(c), 'ME'), "[",2)
colnames(c) = paste0(colnames(c), "_C")


cor_F = t(moduleTraitCor_f)[rows,]
pval_F = t(moduleTraitPvalue_f)[rows,]

f= cor_F
for(i in 1:length(f)){
  f[i] = paste0(cor_F[i], " (", pval_F[i], ")")
}  

colnames(f) = sapply(strsplit(colnames(f), 'ME'), "[",2)
colnames(f) = paste0(colnames(f), "_F")



cor_M = t(moduleTraitCor_m)[rows,]
pval_M = t(moduleTraitPvalue_m)[rows,]

m= cor_M
for(i in 1:length(m)){
  m[i] = paste0(cor_M[i], " (", pval_M[i], ")")
}  

colnames(m) = sapply(strsplit(colnames(m), 'ME'), "[",2)
colnames(m) = paste0(colnames(m), "_M")

S8 = cbind(c,f,m)

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


gtf = rtracklayer::import("results/flat/RNA-seq/mus_stringtie_merged.gtf")
gtf = as.data.frame(gtf)
gtf = gtf[,c("gene_name","ref_gene_id")]
gtf = as.data.frame(unique(gtf))
gtf = gtf[-which(is.na(gtf$ref_gene_id)),]
require(dplyr)
gtf.merged = gtf %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ref_gene_id = paste(ref_gene_id, collapse = ","))

local_eqtl$name2 = sapply(strsplit(local_eqtl$Gene.Name, split = "_isoform"),"[",1)

local_eqtl = merge(local_eqtl, gtf.merged, by.x = "name2", by.y = "gene_name", all.x=T)

local_eqtl = local_eqtl[-c(1,2)]
local_eqtl = local_eqtl[,c(1,12, 2:11)]
colnames(local_eqtl)[2] = "Ensembl_ID"



write.csv(local_eqtl, file = "~/Desktop/supp_tables/S11.csv", row.names = F)

#S12
#qsox1 cluster
library(Seurat)
load("./results/Rdata/seurat_ob.Rdata") #loads as "ob



ob.markers <- FindAllMarkers(ob, only.pos = TRUE)


ob.markers[which(ob.markers$gene =="Qsox1"),] # 11

# cluster 1 genes
cluster1.markers <- FindMarkers(ob, ident.1 = 1, min.pct = 0.25,only.pos = T)
S12 = cluster1.markers[which(cluster1.markers$p_val_adj <= 0.05),]

mart = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations

out <- (getBM(attributes=c('mgi_id','external_gene_name'),
              filters = 'external_gene_name', values = rownames(S12), mart = mart))
S12$gene = rownames(S12)
S12 = merge(S12, out, by.x="gene", by.y = "external_gene_name", all.x=T)

S12[which(is.na(S12$mgi_id)),"mgi_id"] = c("MGI:1919311", "MGI:3642298", "MGI:3649931", "MGI:3606192", "MGI:1328326")

S12 = S12[,c(1,7,2:6)]
colnames(S12)[1:2] = c("Gene_Name", "MGI_ID")
write.csv(S12, file = "~/Desktop/supp_tables/S12.csv", row.names = T)
















#S13 nad S14 from Larry