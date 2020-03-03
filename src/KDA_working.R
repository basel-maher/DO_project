library("igraph")
library("bnlearn")
library("parallel")
library("dplyr")
library(ggplot2)
library(ggpubr)
options(stringsAsFactors = FALSE)

bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}


annot_file = read.csv("~/Documents/projects/DO_project/results/flat/annot_file.csv")
annot_file = annot_file[,c(1,2)]

#load("./results/Rdata/networks/moduleTraitPvalue_f.RData")
#moduleTraitPvalue = moduleTraitPvalue_f

#load("./results/Rdata/networks/moduleTraitPvalue_m.RData")
#moduleTraitPvalue = moduleTraitPvalue_m

#load("./results/Rdata/networks/moduleTraitPvalue_f_BICOR.RData")
#moduleTraitPvalue = moduleTraitPvalue_f


#load("./results/Rdata/networks/moduleTraitPvalue_f_BICOR.RData")
#load("./results/Rdata/networks/moduleTraitPvalue_full_4.RData")

load("./results/Rdata/networks/moduleTraitPvalue_full_4.RData")
#load("./results/Rdata/networks/geneModMemAnnot_f_power8_BICOR.RData")
#load("./results/Rdata/networks/geneModMemAnnot_m_power8_BICOR.RData")
#load("./results/Rdata/networks/geneModMemAnnot_power7.RData")
#load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")
#load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
#load("./results/Rdata/networks/geneModMemAnnot_sexcombined_power7_BICOR.RData")
#use this for sex specific networks, but not full nets
#geneModMemAnnot = combat_annot_f

geneModMemAnnot = combat_annot_m

#gene superset
#superset = read.delim("./results/flat/superduperset_GO_MGI_IMPC.txt", stringsAsFactors = FALSE, header = FALSE)
superset = read.delim("./results/flat/superduperset_sansGWAS.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]

#BNs learned on high performance computing cluster
bns = list.files("~/Documents/projects/DO_project/results/Rdata/networks/bn_4/")
#bns = "hybrid_brown_nobl_4.Rdata"
#
# x = bn2igraph(brown_bn)
# subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,3,nodes = "Nab2"))))
# plot(subgraph,vertex.label.cex=0.65,edge.width=2, vertex.size=20, margin=-0.4, vertex.label.dist=0.2, vertex.label.degree=-pi)

####

out = list()
counter = 1
for(net in bns){
  print(counter)
  color = strsplit(net,"_")[[1]][2]
  #print(color)
  load(paste0("~/Documents/projects/DO_project/results/Rdata/networks/bn_4/",net))
  obj_name = paste0(color,"_bn")
  assign(x = obj_name ,bn)
  
  z = bn2igraph(get(obj_name)) #get "gets" an object from env based on name. otherwise its just a string here
  
  mod_genes = geneModMemAnnot[which(geneModMemAnnot$color == color),"gene"]
  
  if(length(which(is.na(mod_genes))) > 0){
    mod_genes = mod_genes[-which(is.na(mod_genes))]
  }
  
  kda = as.data.frame(matrix(nrow=length(mod_genes),ncol = 8))
  
  
  
  
  i=0
  for(gene in mod_genes){
    
    mod_genes_cur = mod_genes[-which(mod_genes == gene)]
    
    neighbors = names(unlist(neighborhood(z,nodes = gene,order = 3)))
    degree = length(names(unlist(neighborhood(z,nodes = gene,order = 1))))-1
    
    num_neib = length(neighbors)
    
    num_bone_neib = length(which(tolower(neighbors) %in% tolower(superset)))
    in_neib_not_bone = num_neib - num_bone_neib
    not_neib_in_bone = length(which(tolower(mod_genes_cur) %in% tolower(superset) & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
    
    not_neib_not_bone = length(which(tolower(mod_genes_cur) %in% tolower(superset) == FALSE & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
    
    num_in = length(neighbors(z, gene, mode="in")) #0 "in" arcs signifies a root node
    num_out = length(neighbors(z, gene, mode="out"))
    
    i=i+1
    kda[i,1] = gene
    kda[i,2] = color
    kda[i,3] = num_neib
    kda[i,4] = num_bone_neib
    kda[i,5] = num_in
    kda[i,6] = degree
    kda[i,7] = kda[i,4]/kda[i,3]
    kda[i,8] = length(mod_genes_cur)
    kda[i,9] = not_neib_in_bone
    kda[i,10] = in_neib_not_bone
    kda[i,11] = not_neib_not_bone
    
    num_bone_genes_inMod = length(which(tolower(mod_genes_cur) %in% tolower(superset)))
    kda[i,12] = num_bone_genes_inMod
    kda[i,13] = num_out
    #kda[i,11] = x
  }
  
  colnames(kda) = c("gene","color","num_neib","num_bone_neib","num_in_arcs","degree","ratio_bone_neib","num_genes_inMod", "not_neib_in_bone","in_neib_not_bone","not_neib_not_bone","num_bone_genes_inMod", "num_out")
  
  #kda$fdr_pval = p.adjust(kda$phyper,method = "BH")
  
  out[[counter]] = kda
  names(out)[counter] = color
  counter = counter+1
}


kda_full = out



####zhang####
#network size mu
#candidate driver: neib > mu-bar + sd(mu)
for(i in 1:length(kda_full)){#for module in full net
  
  #candidate driver: neib > mu-bar + sd(mu)
  #num neib is 2 step? look at 2 step or degree?
  kda_full[[i]]$driver = 0
  kda_full[[i]]$driver[which(kda_full[[i]]$num_neib > mean(kda_full[[i]]$num_neib) + sd(kda_full[[i]]$num_neib))] = 1
  
  #2 sd
  kda_full[[i]]$driver2 = 0
  kda_full[[i]]$driver2[which(kda_full[[i]]$num_neib > mean(kda_full[[i]]$num_neib) + 2*(sd(kda_full[[i]]$num_neib)))] = 1
  
  #hub genes (d-bar + 2sd(d))
  #num out arcs or degree total?
  kda_full[[i]]$hub = 0
  kda_full[[i]]$hub[which(kda_full[[i]]$num_out > (mean(kda_full[[i]]$num_out) + 2*(sd(kda_full[[i]]$num_out))))] = 1
}
zhang = bind_rows(kda_full)

#zhang = zhang[which(zhang$driver==1 | zhang$hub==1),]
#############
kda_full = zhang

all = bind_rows(kda_full)
all = all[-which(all$num_neib<=2),] #remove unconnected genes or those connected to only 1 

thresh = mean(all$num_neib) - sd(all$num_neib)
all = all[-which(all$num_neib < thresh)]


for(i in 1:nrow(all)){
  all$hyperZhang[i] = phyper(q=all$num_bone_neib[i]-1, m=length(superset), n = nrow(zhang) - length(superset), k=all$num_neib[i], lower.tail = FALSE)
  
  #mod = subset(all, all$color == all$color[i])
  #superMod = length(which(tolower(superset) %in% tolower(mod$gene)))
  #print(superMod)
  #all$hyper_mod[i] = phyper(q=all$num_bone_neib[i]-1, m=superMod, n = nrow(mod) - superMod, k=all$num_neib[i], lower.tail = FALSE)

}

all$hyper_bonf = NA
#all$hyper_mod_bonf = NA
# for(i in which(all$gene %in% zhang$gene)){
#   gene = all$gene[i]
#   zhang[which(zhang$gene == gene),"hyper"] = all$hyperZhang[i]
#   zhang[which(zhang$gene == gene),"hyper_bonf"] = all$hyperZhang_fdr[i]
# }
all$hyper_bonf = p.adjust(all$hyperZhang, method="fdr")
#all$hyper_mod_bonf = p.adjust(all$hyper_mod, method="fdr")
all_sexcombined_power7_BICOR = all



zhang_combined_BICOR = zhang

zhang_males_allModules_MGI_GO_ONLY = zhang


#sig_mod = moduleTraitPvalue[which(rownames(moduleTraitPvalue) %in% names(which(apply(moduleTraitPvalue, 1, function(r) any(r < 0.05/36))))),]
sig_mod = moduleTraitPvalue[which(rownames(moduleTraitPvalue) %in% names(which(apply(moduleTraitPvalue, 1, function(r) any(r < 0.05/length(unique(zhang$color))))))),]
sig_mod = rownames(sig_mod)
sig_mod = unlist(strsplit(sig_mod, "ME"))[seq(2,length(sig_mod)*2,by = 2)]

zhang = zhang[which(zhang$color %in% sig_mod),]
zhang_m_sigModules = zhang

