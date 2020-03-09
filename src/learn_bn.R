library(igraph)
library(bnlearn)
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
#module membership for all genes in all modules (from auto_wgcna.R)
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
