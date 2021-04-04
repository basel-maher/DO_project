library(rtracklayer)
##supp tables

#S1
#from est_herit

#S2 in rev1
####NEW SUPP TABLE WITH ALL MOUSE PHENOTYPIC DATA - RAW
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")
x = colnames(pheno_combined) #from S2 below (Holm-Bonf cor table)

pheno_combined = as.data.frame(cross_basic$pheno)
rownames(pheno_combined) = pheno_combined$Mouse.ID

pheno_combined = pheno_combined[,x]
pheno_combined$Mouse.ID = rownames(pheno_combined)
pheno_combined = pheno_combined[,c(56,1:55)] # mouse ID first

covar = cross_basic$covar
covar$Mouse.ID = rownames(covar)

S0 = merge(covar, pheno_combined, by="Mouse.ID")
S0 = S0[order(as.numeric(S0$Mouse.ID)),]

S0 = S0[,-c(8:17)]


write.csv(S0, file = "~/Desktop/nat_com_revs/supp/S2_rev.csv")



#S3
#Holm-Bonferroni correlations with femoral strength
library(psych)
set.seed(8675309)
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
#x = psych::corr.test(pheno_combined,method = "p")
xs = psych::corr.test(pheno_combined, method = "s", use="p")
#p vals
#p = x$p[,"bending_max_load"]
ps = signif(xs$p[,"bending_max_load"],3)
#correlation
#r = x$r[,"bending_max_load"]
rs = signif(xs$r[,"bending_max_load"],3)

#S.E.
#se = x$se[,"bending_max_load"]
ses = signif(xs$se[,"bending_max_load"],3)

table = cbind(rs,ps,ses)


write.csv(table, file = "~/Desktop/nat_com_revs/supp/S3_rev.csv")



#S4
#table = x$r
#p = x$p

#for( i in 1:length(table)){
#  table[i] = paste0(signif(as.numeric(table[i],3)), "(", signif(as.numeric(p[i],3)), ")")
#}

table2 = xs$r
p = xs$p

for( i in 1:length(table2)){
  table2[i] = paste0(signif(as.numeric(table2[i]),3), " (", signif(as.numeric(p[i]),3), ")")
}

write.csv(table2, file = "~/Desktop/nat_com_revs/supp/S4_rev.csv")


#S5

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

require(data.table)
MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)
S4$name2 = tolower(sapply(strsplit(S4$Gene.Name, split = "_isoform"),"[",1))

S4 = merge(S4, MGI_markers, by.x="name2", by.y="name", all.x=T)

idx = which(is.na(S4$`MGI Accession ID`))

for(i in idx){
  
  zz = grep(x = MGI_markers$`Marker Synonyms (pipe-separated)`, pattern = paste0("\\b",S4$name2[i],"\\b"),ignore.case = T)
  
  if(length(zz) == 1){
    S4$`MGI Accession ID`[i] = MGI_markers$`MGI Accession ID`[zz]
    print(paste0(S4$name2[i], " ----- ", MGI_markers$`Marker Synonyms (pipe-separated)`[zz]))
  }
  
}
S4[which(is.na(S4$`MGI Accession ID`)),"name2"]
S4[which(S4$name2 == "grasp"),"MGI Accession ID"] = "MGI:1860303"
S4[which(S4$name2 == "marc2"),"MGI Accession ID"] = "MGI:1914497"

gtf = rtracklayer::import("results/flat/RNA-seq/mus_stringtie_merged.gtf")
gtf = as.data.frame(gtf)
gtf = gtf[,c("gene_name","ref_gene_id")]
gtf = as.data.frame(unique(gtf))
gtf = gtf[-which(is.na(gtf$ref_gene_id)),]
require(dplyr)
gtf.merged = gtf %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ref_gene_id = paste(ref_gene_id, collapse = ","))

S4$alt_id = NA
idx = which(is.na(S4$`MGI Accession ID`))

for(i in idx){
  
  zz = which(tolower(gtf.merged$gene_name) == S4$name2[i])
  
  if(length(zz) == 1){
    S4$alt_id[i] = gtf.merged$ref_gene_id[zz]
    print(paste0(S4$name2[i], " ----- ", gtf.merged$gene_name[zz]))
  }
  
}

S4 = S4[,c(2,130,142,3:129)]

S4[,7:ncol(S4)] = signif(S4[,7:ncol(S4)],3)

write.csv(S4, file = "~/Desktop/nat_com_revs/supp/S5_rev.csv", row.names = F)




#S6 - top GO terms for each network, sorted by pval
library(GO.db)
goterms = Term(GOTERM)

load("./results/Rdata/networks/GO_sft4.RData")

networks = as.data.frame(matrix(nrow=10000, ncol = 10000))
nrow_counter = c()
net_counter = 1
col_counter = 1
for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  n = n[which(n$classic <=0.05),]
  
  x = n$GO.ID
  net_terms = unname(goterms[match(x, names(goterms))])
  
  networks[1:nrow(n),col_counter] = net_terms
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_C_GO_term")
  
  col_counter = col_counter + 1

  networks[1:nrow(n),col_counter] = n$classic
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_C_p_value")
  nrow_counter = append(nrow_counter, nrow(n))
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}




load("./results/Rdata/networks/GO_Females_sft4.RData")
net_counter = 1

for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  n = n[which(n$classic <=0.05),]
  
  x = n$GO.ID
  net_terms = unname(goterms[match(x, names(goterms))])
  
  
  networks[1:nrow(n),col_counter] = net_terms
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_F_GO_term")
  
  col_counter = col_counter + 1
  
  networks[1:nrow(n),col_counter] = n$classic
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_F_p_value")
  nrow_counter = append(nrow_counter, nrow(n))
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}





load("./results/Rdata/networks/GO_Males_sft5.RData")


net_counter = 1

for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  n = n[which(n$classic <=0.05),]
  
  x = n$GO.ID
  net_terms = unname(goterms[match(x, names(goterms))])
  
  
  networks[1:nrow(n),col_counter] = net_terms
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_M_GO_term")
  
  col_counter = col_counter + 1
  
  networks[1:nrow(n),col_counter] = n$classic
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_M_p_value")
  nrow_counter = append(nrow_counter, nrow(n))
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}


networks = networks[c(1:max(nrow_counter)),c(1:248)]



write.csv(networks, file = "~/Desktop/nat_com_revs/supp/S6_rev.csv", row.names = F,na = "-")


#S7
#these are mouse genes
superset = read.delim("./results/flat/superset_in_networks.txt", stringsAsFactors = FALSE, header = FALSE)

#superset = superset[,1]
#superset = capitalize(superset)

require(data.table)
MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)

superset = merge(superset, MGI_markers, by.x="V1", by.y="name", all.x=T)

#which(duplicated(superset$V1)) ##no dups anymore
#superset = superset[-1356,]

zz = which(is.na(superset$`MGI Accession ID`) == FALSE)

superset$MGI_NAME[zz] = superset$`Marker Symbol`[zz]

zz = which(is.na(superset$`MGI Accession ID`))
require(Hmisc)
superset$MGI_NAME[zz] = capitalize(superset$V1[zz])

#superset[which(is.na(superset$`MGI Accession ID`)),"MGI Accession ID"] = c("MGI:2140364","MGI:3575190", "MGI:3575247", "MGI:3575248", "MGI:5573174","MGI:1915720","MGI:3819962", "MGI:5633762","MGI:3812132")

#superset[which(superset$MGI_NAME == "Adprhl2"),"Feature Type"] = "protein coding gene"
superset[which(superset$MGI_NAME == "Adprhl2"),] = c("adprhl2","MGI:2140364","4","60.36","126316351","126321703","-","Adprs","O", "ADP-ribosylserine hydrolase","Gene","protein coding gene","Adprhl2 Arh3","Adprhl2")
#superset[which(superset$MGI_NAME == "Impad1"),"Feature Type"] = "protein coding gene"
superset[which(superset$MGI_NAME == "Impad1"),] = c("impad1","MGI:1915720","4","2.56","4762484","4793306","-","Bpnt2","O","3'(2'), 5'-bisphosphate nucleotidase 2","Gene","protein coding gene","gPAPP 1110001C20Rik Jaws Impad1","Impad1")
#superset[which(superset$MGI_NAME %in% c("Hlb324b",  "Hlb328",   "Hlb330",   "Hydro",   "M1665asr", "Rdns","Shsn")),"Feature Type"] = "heritable phenotypic marker"

#superset_pruned = superset[which(superset$`Feature Type` %in% c("complex/cluster/region","heritable phenotypic marker","QTL","unclassified other genome feature")==FALSE),14]
superset=superset[,c(14,2)]

write.csv(superset, file = "~/Desktop/nat_com_revs/supp/S7_rev.csv", row.names = F)




# full_net = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv", stringsAsFactors = FALSE)
# female_net = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv", stringsAsFactors = FALSE)
# male_net = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv", stringsAsFactors = FALSE)
# 
# merge(full_net, female_net, male_net, )


#S8
full = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4_REV.csv", stringsAsFactors = F)
full = full[,c(1:4,16,17)]

colnames(full)[2:ncol(full)] = paste0(colnames(full)[2:ncol(full)], "_C")



female = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4_REV.csv", stringsAsFactors = F)
female = female[,c(1:4,17,18)]

colnames(female)[2:ncol(female)] = paste0(colnames(female)[2:ncol(female)], "_F")


male = read.csv("./results/flat/key_driver_analysis_MALES_sft5_REV.csv", stringsAsFactors = F)
male = male[,c(1:4,17,18)]

colnames(male)[2:ncol(male)] = paste0(colnames(male)[2:ncol(male)], "_M")

S7 = merge(full, female, by="gene", all=T)
S7 = merge(S7, male, by="gene", all=T)

S7 = S7[,c(1,2,7,12,3:6,8:11,13:16)]
colnames(S7) = c("gene","module_C","module_F","module_M","num_neib_C","num_bone_neib_C","nominal_pval_C", "FDR_pval_C", "num_neib_F","num_bone_neib_F","nominal_pval_F","FDR_pval_F","num_neib_M","num_bone_neib_M","nominal_pval_M","FDR_pval_M")

require(data.table)
MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)
S7$name2 = tolower(sapply(strsplit(S7$gene, split = "_isoform"),"[",1))

S7 = merge(S7, MGI_markers, by.x="name2", by.y="name", all.x=T)

idx = which(is.na(S7$`MGI Accession ID`))

for(i in idx){
  
  zz = grep(x = MGI_markers$`Marker Synonyms (pipe-separated)`, pattern = paste0("\\b",S7$name2[i],"\\b"),ignore.case = T)
  
  if(length(zz) == 1){
    S7$`MGI Accession ID`[i] = MGI_markers$`MGI Accession ID`[zz]
    print(paste0(S7$name2[i], " ----- ", MGI_markers$`Marker Synonyms (pipe-separated)`[zz]))
  }
  
}

S7[which(is.na(S7$`MGI Accession ID`)),"name2"]
S7[which(S7$name2 == "grasp"),"MGI Accession ID"] = "MGI:1860303"
S7[which(S7$name2 == "marc2"),"MGI Accession ID"] = "MGI:1914497"

gtf = rtracklayer::import("results/flat/RNA-seq/mus_stringtie_merged.gtf")
gtf = as.data.frame(gtf)
gtf = gtf[,c("gene_name","ref_gene_id")]
gtf = as.data.frame(unique(gtf))
gtf = gtf[-which(is.na(gtf$ref_gene_id)),]
require(dplyr)
gtf.merged = gtf %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ref_gene_id = paste(ref_gene_id, collapse = ","))

S7$alt_id = NA
idx = which(is.na(S7$`MGI Accession ID`))

for(i in idx){
  
  zz = which(tolower(gtf.merged$gene_name) == S7$name2[i])
  
  if(length(zz) == 1){
    S7$alt_id[i] = gtf.merged$ref_gene_id[zz]
    print(paste0(S7$name2[i], " ----- ", gtf.merged$gene_name[zz]))
  }
  
}


S7 = S7[,c(2,18,30,3:17)]


#get BANs and add BAN column
all = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4_REV.csv", stringsAsFactors = F)
all_f = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4_REV.csv", stringsAsFactors = F)
all_m = read.csv("./results/flat/key_driver_analysis_MALES_sft5_REV.csv", stringsAsFactors = F)

BANs = c(all[which(all$hyper<=0.05), "gene"], all_m[which(all_m$hyper<=0.05), "gene"], all_f[which(all_f$hyper<=0.05), "gene"])
BANs = unique(BANs)



idx = which(S7$gene %in% BANs)
S7$BAN = "NO"
S7$BAN[idx] = "YES"

S7[,c("nominal_pval_C","FDR_pval_C","nominal_pval_F","FDR_pval_F","nominal_pval_M","FDR_pval_M")] = signif(S7[,c("nominal_pval_C","FDR_pval_C","nominal_pval_F","FDR_pval_F","nominal_pval_M","FDR_pval_M")],3)

write.csv(S7, file = "~/Desktop/nat_com_revs/supp/S8_rev.csv", row.names = F)


#S9
fn = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_FNBMD_over75_REV.txt")
ls = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_LSBMD_over75_REV.txt")
morris_coloc = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_morris_v7_all_results_over75_REV.txt")

S9 = rbind(morris_coloc, fn, ls)
S9 = S9[which(S9$H4 >= 0.75),]

S9[which(S9$pheno == "BMD"),"pheno"] = "eBMD"

S9[,c("H0","H1","H2","H3","H4")] = signif(S9[,c("H0","H1","H2","H3","H4")],3)

write.csv(S9, file = "~/Desktop/nat_com_revs/supp/S9_rev.csv", row.names = F)


#S10 correlation between modules and bone phenotypes
load(file = "./results/Rdata/networks/moduleTraitPvalue_full_4.RData")
load(file = "./results/Rdata/networks/moduleTraitCor_full_4.RData")

load(file = "./results/Rdata/networks/moduleTraitPvalue_m.RData")
load(file = "./results/Rdata/networks/moduleTraitCor_m.RData")

load(file = "./results/Rdata/networks/moduleTraitPvalue_f.RData")
load(file = "./results/Rdata/networks/moduleTraitCor_f.RData")

rows = match(colnames(pheno_combined), colnames(moduleTraitCor))

cor = t(moduleTraitCor)[rows,]
cor = signif(cor,3)
pval = t(moduleTraitPvalue)[rows,]
pval = signif(pval,3)

c = cor
for(i in 1:length(c)){
  c[i] = paste0(cor[i], " (", pval[i], ")")
}  

colnames(c) = sapply(strsplit(colnames(c), 'ME'), "[",2)
colnames(c) = paste0(colnames(c), "_C")


cor_F = t(moduleTraitCor_f)[rows,]
cor_F = signif(cor_F,3)
pval_F = t(moduleTraitPvalue_f)[rows,]
pval_F = signif(pval_F,3)

f= cor_F
for(i in 1:length(f)){
  f[i] = paste0(cor_F[i], " (", pval_F[i], ")")
}  

colnames(f) = sapply(strsplit(colnames(f), 'ME'), "[",2)
colnames(f) = paste0(colnames(f), "_F")



cor_M = t(moduleTraitCor_m)[rows,]
cor_M = signif(cor_M,3)
pval_M = t(moduleTraitPvalue_m)[rows,]
pval_M = signif(pval_M,3)

m= cor_M
for(i in 1:length(m)){
  m[i] = paste0(cor_M[i], " (", pval_M[i], ")")
}  

colnames(m) = sapply(strsplit(colnames(m), 'ME'), "[",2)
colnames(m) = paste0(colnames(m), "_M")

S8 = cbind(c,f,m)



write.csv(S8, file = "~/Desktop/nat_com_revs/supp/S10_rev.csv", row.names = T)


#11

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
write.csv(S9, file = "~/Desktop/nat_com_revs/supp/S11_rev.csv", row.names = F)






#S12
#SIFT annotations
sift = read.table("~/Desktop/nat_com_revs/supp/VEP_output_050320.txt", header=T)
write.csv(sift, file = "~/Desktop/nat_com_revs/supp/S12_rev.csv", row.names = F, quote=T)


#S13
#all local eqtl

load("./results/Rdata/local_eqtl.Rdata")
local_eqtl = local_eqtl[,c(1,8,3,4,5,6,7,9,10,11,12,16)]

head(local_eqtl)
local_eqtl$Start = local_eqtl$Start/1000000
local_eqtl$End = local_eqtl$End/1000000


require(data.table)
MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)
local_eqtl$name2 = tolower(sapply(strsplit(local_eqtl$Gene.Name, split = "_isoform"),"[",1))

local_eqtl = merge(local_eqtl, MGI_markers, by.x="name2", by.y="name", all.x=T)

idx = which(is.na(local_eqtl$`MGI Accession ID`))

for(i in idx){
  
  zz = grep(x = MGI_markers$`Marker Synonyms (pipe-separated)`, pattern = paste0("\\b",local_eqtl$name2[i],"\\b"),ignore.case = T)
  
  if(length(zz) == 1){
    local_eqtl$`MGI Accession ID`[i] = MGI_markers$`MGI Accession ID`[zz]
    print(paste0(local_eqtl$name2[i], " ----- ", MGI_markers$`Marker Synonyms (pipe-separated)`[zz]))
  }
  
}

local_eqtl[which(is.na(local_eqtl$`MGI Accession ID`)),"name2"]
local_eqtl[which(local_eqtl$name2 == "grasp"),"MGI Accession ID"] = "MGI:1860303"
local_eqtl[which(local_eqtl$name2 == "marc2"),"MGI Accession ID"] = "MGI:1914497"
local_eqtl[which(local_eqtl$name2 == "d930016d06rik"),"MGI Accession ID"] = "MGI:3040669"

gtf = rtracklayer::import("results/flat/RNA-seq/mus_stringtie_merged.gtf")
gtf = as.data.frame(gtf)
gtf = gtf[,c("gene_name","ref_gene_id")]
gtf = as.data.frame(unique(gtf))
gtf = gtf[-which(is.na(gtf$ref_gene_id)),]
require(dplyr)
gtf.merged = gtf %>%
  dplyr::group_by(gene_name) %>%
  dplyr::summarise(ref_gene_id = paste(ref_gene_id, collapse = ","))

local_eqtl$alt_id = NA
idx = which(is.na(local_eqtl$`MGI Accession ID`))

for(i in idx){
  
  zz = which(tolower(gtf.merged$gene_name) == local_eqtl$name2[i])
  
  if(length(zz) == 1){
    local_eqtl$alt_id[i] = gtf.merged$ref_gene_id[zz]
    print(paste0(local_eqtl$name2[i], " ----- ", gtf.merged$gene_name[zz]))
  }
  
}

local_eqtl = local_eqtl[,c(3,14,26,4:13)]

write.csv(local_eqtl, file = "~/Desktop/nat_com_revs/supp/S13_rev.csv", row.names = F)



#S14
#qsox1 cluster
library(Seurat)
load("./results/Rdata/seurat_ob.Rdata") #loads as "ob



ob.markers <- FindAllMarkers(ob, only.pos = TRUE)


ob.markers[which(ob.markers$gene =="Qsox1"),] # 11

# cluster 1 genes
cluster1.markers <- FindMarkers(ob, ident.1 = 1, min.pct = 0.25,only.pos = T)
S12 = cluster1.markers[which(cluster1.markers$p_val_adj <= 0.05),]

S12$name = rownames(S12)

require(data.table)
MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)
S12$name2 = tolower(sapply(strsplit(S12$name, split = "_isoform"),"[",1))

S12 = merge(S12, MGI_markers, by.x="name2", by.y="name", all.x=T)

idx = which(is.na(S12$`MGI Accession ID`))

for(i in idx){
  
  zz = grep(x = MGI_markers$`Marker Synonyms (pipe-separated)`, pattern = paste0("\\b",S12$name2[i],"\\b"),ignore.case = T)
  
  if(length(zz) == 1){
    S12$`MGI Accession ID`[i] = MGI_markers$`MGI Accession ID`[zz]
    print(paste0(S12$name2[i], " ----- ", MGI_markers$`Marker Synonyms (pipe-separated)`[zz]))
  }
  
}

write.csv(S12, file = "~/Desktop/nat_com_revs/supp/S14_rev.csv", row.names = F)
















#S13 nad S14 from Larry