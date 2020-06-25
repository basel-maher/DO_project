library(igraph)
library(bnlearn)
########
#Run on supercomputaing cluster
########
load("./edata_4.RData")
edata=edata_trim
load("./geneModMemAnnot_power4.RData")

#BASH STUFF#
args <- commandArgs()
print(args)
num = as.numeric(args[6])
print(num)
###

resid = edata
connect = geneModMemAnnot

clr = unique(connect$color)

color_a = clr[num]
print(color_a)

mod_genes = connect[which(connect$color==color_a),"gene"]
mod_genes_membership = connect[which(connect$gene %in% mod_genes),c("gene","gene_id",paste0("ME",color_a),"color")]
mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == color_a),]
mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$gene)])

bn = mmhc(mod_genes_exp)
print("check")
varName = paste0("hybrid_",color_a,"_nobl_4")

save(bn, file = paste0("./bn_4/",varName,".RData"))



########Females

load("./edata_trim_f_4.RData")
resid=edata_trim_f
load("./geneModMemAnnot_f_power4.RData")
connect = combat_annot_f

#BASH STUFF#
args <- commandArgs()
print(args)
num = as.numeric(args[6])
print(num)
###


clr = unique(connect$color)

color_a = clr[num]
print(color_a)

mod_genes = connect[which(connect$color==color_a),"gene"]
mod_genes_membership = connect[which(connect$gene %in% mod_genes),c("gene","Gene.Name",paste0("ME",color_a),"color")]
mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == color_a),]

mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$Gene.Name)])


bn = mmhc(mod_genes_exp)
print("check")
varName = paste0("hybrid_",color_a,"_f_4")

save(bn, file = paste0("./bn_f_4/",varName,".RData"))



#####Males
load("./edata_trim_m_5.RData")
resid=edata_trim_m
load("./geneModMemAnnot_m_power5.RData")

connect = combat_annot_m


#BASH STUFF#
args <- commandArgs()
print(args)
num = as.numeric(args[6])
print(num)
###


clr = unique(connect$color)

color_a = clr[num]
print(color_a)

mod_genes = connect[which(connect$color==color_a),"gene"]
mod_genes_membership = connect[which(connect$gene %in% mod_genes),c("gene","Gene.Name",paste0("ME",color_a),"color")]
mod_genes_membership = mod_genes_membership[which(mod_genes_membership$color == color_a),]

mod_genes_exp = as.data.frame(resid[,which(colnames(resid) %in% mod_genes_membership$Gene.Name)])



bn = mmhc(mod_genes_exp)
print("check")
varName = paste0("hybrid_",color_a,"_m_5")

save(bn, file = paste0("./bn_m_5/",varName,".RData"))

