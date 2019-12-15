library("igraph")
library("bnlearn")
library("parallel")
library("dplyr")
library(ggplot2)
library(ggpubr)
options(stringsAsFactors = FALSE)

#load("./results/Rdata/networks/moduleTraitPvalue_f.RData")
#moduleTraitPvalue = moduleTraitPvalue_f

#load("./results/Rdata/networks/moduleTraitPvalue_f.RData")
#moduleTraitPvalue = moduleTraitPvalue_f
#load("./results/Rdata/networks/moduleTraitPvalue_f.RData")



load("./results/Rdata/networks/moduleTraitPvalue_full_4.RData")

#load("./results/Rdata/networks/geneModMemAnnot_power7.RData")
#load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")
#load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
#load("./results/Rdata/networks/geneModMemAnnot_power7_BICOR.RData")
#use this for sex specific networks, but not full nets
#geneModMemAnnot = combat_annot_f

geneModMemAnnot = combat_annot_f

#required to make plots and such. conver bn to igraph object
bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}

#gene superset
#superset = read.delim("./results/flat/superduperset_GO_MGI_IMPC.txt", stringsAsFactors = FALSE, header = FALSE)
superset = read.delim("./results/flat/superduperset_sansGWAS.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]

#BNs learned on high performance computing cluster
bns = list.files("~/Desktop/bn_4/")
#bns = "hybrid_brown_nobl_4.Rdata"
#


x = bn2igraph(brown_bn)
subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,3,nodes = "Nab2"))))
plot(subgraph,vertex.label.cex=0.65,edge.width=2, vertex.size=20, margin=-0.4, vertex.label.dist=0.2, vertex.label.degree=-pi)

####

out = list()
counter = 1
for(net in bns){
  print(counter)
  color = strsplit(net,"_")[[1]][2]
  #print(color)
  load(paste0("~/Desktop/bn_4/",net))
  obj_name = paste0(color,"_bn")
  assign(x = obj_name ,bn)
  
  z = bn2igraph(get(obj_name)) #get "gets" an object from env based on name. otherwise its just a string here
  mod_genes = geneModMemAnnot[which(geneModMemAnnot$color == color),"gene"]
  
  if(length(which(is.na(mod_genes))) > 0){
    mod_genes = mod_genes[-which(is.na(mod_genes))]
  }
  
  #mod_genes = colnames(resid)[which(colnames(resid) %in% mod_genes)]
  ##when to drop lowly connected nodes?
  #low = arcs(get(obj_name))[which(neighborhood.size(z,nodes = arcs(get(obj_name)),order = 3)==2)]
  #mod_genes = mod_genes[-which(mod_genes %in% low)]

  kda = as.data.frame(matrix(nrow=length(mod_genes),ncol = 8))
  
  i=0
  for(gene in mod_genes){
    
    mod_genes_cur = mod_genes[-which(mod_genes == gene)]
    
    neighbors = names(unlist(neighborhood(z,nodes = gene,order = 2)))
    degree = length(names(unlist(neighborhood(z,nodes = gene,order = 1))))-1
    
    #neighbors = neighbors[-which(neighbors == gene)]
    
    num_neib = length(neighbors)
    
    num_bone_neib = length(which(tolower(neighbors) %in% tolower(superset)))
    in_neib_not_bone = num_neib - num_bone_neib
    not_neib_in_bone = length(which(tolower(mod_genes_cur) %in% tolower(superset) & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))

    not_neib_not_bone = length(which(tolower(mod_genes_cur) %in% tolower(superset) == FALSE & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
    
    num_in = length(neighbors(z, gene, mode="in")) #0 "in" arcs signifies a root node
    
    i=i+1
    kda[i,1] = gene
    kda[i,2] = color
    kda[i,3] = num_neib
    kda[i,4] = num_bone_neib
    kda[i,5] = num_in
    kda[i,6] = degree
    kda[,7] = kda[,4]/kda[,3]
    kda[,8] = length(mod_genes_cur)
    kda[,9] = not_neib_in_bone
    kda[,10] = in_neib_not_bone
    kda[,11] = not_neib_not_bone
    
    num_bone_genes_inMod = length(which(tolower(mod_genes_cur) %in% tolower(superset)))
    kda[,12] = num_bone_genes_inMod
    
    #kda[i,11] = x
  }
  
  colnames(kda) = c("gene","color","num_neib","num_bone_neib","num_in_arcs","degree","ratio_bone_neib","num_genes_inMod", "not_neib_in_bone","in_neib_not_bone","not_neib_not_bone","num_bone_genes_inMod")
  
  #kda$fdr_pval = p.adjust(kda$phyper,method = "BH")
  
  out[[counter]] = kda
  names(out)[counter] = color
  counter = counter+1
}
#kda$phyper = phyper(q=kda$num_bone_neib-1, m=kda$num_bone_genes, n = kda$num_genes_inMod - kda$num_bone_genes_inMod, k=kda$num_neib, lower.tail = FALSE)
#kda = kda[-which(kda$num_bone_neib==0),]
#kda$phyper_adj = p.adjust(kda$phyper, method = "bonferroni")

kda_full = out



####zhang####
#network size mu
#candidate driver: neib > mu-bar + sd(mu)
for(i in 1:length(kda_full)){
  
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
  kda_full[[i]]$hub[which(kda_full[[i]]$num_out_arcs > mean(kda_full[[i]]$num_out_arcs) + 2*(sd(kda_full[[i]]$num_out_arcs)))] = 1
}
zhang = bind_rows(kda_full)

#zhang = zhang[which(zhang$driver==1 | zhang$hub==1),]
#############
kda_full = out

for(i in 1:length(kda_full)){
  kda_full[[i]] = kda_full[[i]][-which(kda_full[[i]]$num_bone_neib == 0),]
  thresh = mean(kda_full[[i]]$num_neib) - sd(kda_full[[i]]$num_neib)
  kda_full[[i]] = kda_full[[i]][-which(kda_full[[i]]$num_neib < thresh),]
  kda_full[[i]]$hyper = phyper(q=kda_full[[i]]$num_bone_neib, m=kda_full[[i]]$num_bone_genes, n = kda_full[[i]]$num_genes_inMod -  kda_full[[i]]$num_bone_genes_inMod, k=kda_full[[i]]$num_neib, lower.tail = FALSE)
  #kda_full[[i]]$hyperZhang = phyper(q=kda_full[[i]]$num_bone_neib-1, m=1631, n = 16816-1631, k=kda_full[[i]]$num_neib, lower.tail = FALSE)#includes self
  kda_full[[i]]$hyperZhang = phyper(q=kda_full[[i]]$num_bone_neib-1, m=length(superset), n = nrow(zhang)-length(superset), k=kda_full[[i]]$num_neib, lower.tail = FALSE)#includes self
  
  ####SHOULD n BE 16816 - 1631
  #subtract 1 for prob P[X>=x] (null hyp)
  #m=kda_full[[i]]$num_bone_genes
  kda_full[[i]]$hyper_fdr = p.adjust(kda_full[[i]]$hyper, method="bonferroni")
  kda_full[[i]]$hyperZhang_fdr = p.adjust(kda_full[[i]]$hyperZhang, method="bonferroni")
  }

all = bind_rows(kda_full)

#write.csv(all, file="~/Desktop/kda_updatedBoneSet_sansLowlyConn.csv",quote = FALSE)
zhang$hyper = NA
zhang$hyper_bonf = NA

for(i in which(all$gene %in% zhang$gene)){
  gene = all$gene[i]
  zhang[which(zhang$gene == gene),"hyper"] = all$hyperZhang[i]
  zhang[which(zhang$gene == gene),"hyper_bonf"] = all$hyperZhang_fdr[i]
}


##fishers exact test
#contingency table is, for a gene, neighbor and in superset, neighbor not in superset, not neighbor in super, not neighbor not super
# gene = "Cpe" 
# num_neib = neighborhood.size(z,nodes = gene,order = 2)
# num_bone_neib = length(which(tolower(neighbors) %in% tolower(superset[,1])))
# in_neib_not_bone = num_neib - num_bone_neib
# not_neib_in_bone = length(which(tolower(superset[,1]) %in% tolower(neighbors) == FALSE))
# not_neib_not_bone = length(which(tolower(colnames(blue_genes_exp)) %in% tolower(superset[,1]) == FALSE & tolower(colnames(blue_genes_exp)) %in% tolower(neighbors) == FALSE))
# 
# x = matrix(c(num_bone_neib,not_neib_in_bone,in_neib_not_bone,not_neib_not_bone),nrow = 2,ncol = 2)
# fisher.test(x)
###
#how much do we expect for a random gene set? sample from all genes in rna-seq data, same length as num of genes in superset
#two ways. one is sample from all genes in rnaseq data, same length as superset. randomize connections per node.
# second way is to randomly change connections to a gene within a module, using only genes in the module
#annot_file = read.delim("367-FarberDO2_S20.gene_abund.tab",header = TRUE)
#allGenes = unique(annot_file$Gene.Name)
#allGenes = colnames(resid)

#method 2
# out_perm = list()
# counter = 1
# for(net in bns){
#   print(counter)
#   color = strsplit(net,"_")[[1]][2]
#   print(color)
#   load(paste0("~/Desktop/",net))
#   obj_name = paste0(color,"_bn")
#   assign(x = obj_name ,bn)
# 
#   z = bn2igraph(get(obj_name)) #get "gets" an object from env based on name. otherwise its just a string here
#   mod_genes = geneModMemAnnot[which(geneModMemAnnot$color == color),"gene"]
#   perm=0
#   allGenes = geneModMemAnnot$gene
#   m=c()
#   pval_list = list()
# 
#   for(gene in mod_genes){
#     pval_list[[gene]] = vector()
#   }
# 
# 
# 
# 
#   while(perm < 1000){
#     print(perm)
#     kda_rand = as.data.frame(matrix(nrow=length(mod_genes),ncol = 7))
# 
#     colnames(kda_rand) = c("gene","num_neib","num_bone_neib","num_in_arcs","ratio_bone_neib","fisher_pval")
#     kda_rand = as.data.frame(matrix(nrow=length(mod_genes),ncol = 7))
#     i=0
#     for(gene in mod_genes){
#       mod_genes_cur = mod_genes[-which(mod_genes == gene)]
# 
#       num_neib = neighborhood.size(z,nodes = gene,order = 2) - 1
# 
#       neighbors = sample(mod_genes_cur,size = num_neib)
# 
#       num_bone_neib = length(which(tolower(neighbors) %in% tolower(superset)))
#       in_neib_not_bone = num_neib - num_bone_neib
#       not_neib_in_bone = length(which(tolower(mod_genes_cur) %in% tolower(superset) & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
#       not_neib_not_bone = length(which(tolower(mod_genes_cur) %in% tolower(superset) == FALSE & tolower(mod_genes_cur) %in% tolower(neighbors) == FALSE))
#       num_in = length(neighbors(z, gene, mode="in")) #0 "in" arcs signifies a root node
# 
#       i=i+1
#       kda_rand[i,1] = gene
#       kda_rand[i,2] = num_neib
#       kda_rand[i,3] = num_bone_neib
#       kda_rand[i,4] = num_in
#       kda_rand[,5] = kda_rand[,3]/kda_rand[,2]
# 
# 
#       x = matrix(c(num_bone_neib,not_neib_in_bone,in_neib_not_bone,not_neib_not_bone),nrow = 2,ncol = 2)
#       xx = fisher.test(x,alternative = "greater")
#       kda_rand[i,6] = xx$p.value
#       kda_rand[i,7] = kda_rand[i,6]*kda_rand[i,2] # bonf pval = pval * num_neib
#       pval_list[[gene]] =  append(pval_list[[gene]],kda_rand[i,3]) #append bonf (append num bone neib)
# 
#     }
# 
# 
#   perm = perm+1
#   }
# 
#   out_perm[[counter]] = pval_list
#   names(out_perm)[counter] = color
#   counter = counter+1
# }
# 
# save(out_perm, file="~/Desktop/KDA_1000_perms_minus1.Rdata")
# load("~/Desktop/KDA_1000_perms.Rdata")
# 
zhang_combined_BICOR = zhang

zhang_males_allModules_MGI_GO_ONLY = zhang


#sig_mod = moduleTraitPvalue[which(rownames(moduleTraitPvalue) %in% names(which(apply(moduleTraitPvalue, 1, function(r) any(r < 0.05/36))))),]
sig_mod = moduleTraitPvalue[which(rownames(moduleTraitPvalue) %in% names(which(apply(moduleTraitPvalue, 1, function(r) any(r < 0.05/length(unique(zhang$color))))))),]
sig_mod = rownames(sig_mod)
sig_mod = unlist(strsplit(sig_mod, "ME"))[seq(2,length(sig_mod)*2,by = 2)]

zhang = zhang[which(zhang$color %in% sig_mod),]
zhang_m_sigModules = zhang

#######
#some analysis#
#
#are genes that colocalize more highly connected with bone genes in the network?
##basically, do colocalizing genes have more bone neighbors?
coloc_genes = zhang[which(zhang$coloc_H4>=0.75),]
not_coloc_genes = zhang[which(zhang$gene %in% coloc_genes$gene == FALSE),]

wilcox.test(coloc_genes$num_bone_neib, not_coloc_genes$num_bone_neib, alternative = "greater")
#MAYBE THE WAY TO DO THIS IS WITH THE PHYPER RESULTS. DO COLOC GENES HAVE MORE SIGNIFICANT P VALS
wilcox.test(coloc_genes$hyper_bonf, not_coloc_genes$hyper_bonf, alternative = "less")
######
##same but subtract 1 if bone gene because we dont want to count self
coloc_genes$num_bone_neib[which(tolower(coloc_genes$gene) %in% tolower(superset))]  = coloc_genes$num_bone_neib[which(tolower(coloc_genes$gene) %in% tolower(superset))] - 1
not_coloc_genes$num_bone_neib[which(tolower(not_coloc_genes$gene) %in% tolower(superset))]  = not_coloc_genes$num_bone_neib[which(tolower(not_coloc_genes$gene) %in% tolower(superset))] - 1

wilcox.test(coloc_genes$num_bone_neib, not_coloc_genes$num_bone_neib, alternative = "greater")
wilcox.test(coloc_genes$hyper_bonf, not_coloc_genes$hyper_bonf, alternative = "less")

wilcox.test(coloc_genes$ratio_bone_neib, not_coloc_genes$ratio_bone_neib, alternative = "greater")

##are bone genes more highly connected in the network?
x = zhang
x = x[-which(x$num_neib<2),]
bone_df = x[which(tolower(x$gene) %in% tolower(superset)),]
not_bone_df = x[which(tolower(x$gene) %in% tolower(superset) == FALSE),]


wilcox.test(bone_df$degree, not_bone_df$degree, alternative = "greater")#3E-6
wilcox.test(bone_df$num_neib, not_bone_df$num_neib, alternative = "greater")#1E-8
wilcox.test(bone_df$ratio_bone_neib, not_bone_df$ratio_bone_neib, alternative = "greater")#1E-8

x$bone = 0
x[which(tolower(x$gene) %in% tolower(superset)),"bone"] = 1
boxplot(num_neib~bone, data=x)
#

######
#are bone genes more likely to be drivers
bone_df = zhang[which(tolower(zhang$gene) %in% tolower(superset)),]
not_bone_df = zhang[which(tolower(zhang$gene) %in% tolower(superset) == FALSE),]

#bone gene AND drivers / bone genes
#not bone genes AND drivers / not bone genes

num_drivers = c(length(which(bone_df$driver==1)), length(which(not_bone_df$driver == 1)))
num_genes = c(length(bone_df$gene), length(not_bone_df$gene))
prop.test(num_drivers, num_genes, alternative = "greater")
#what if were testing between bone genes and ALL genes
num_drivers = c(length(which(bone_df$driver==1)), length(which(zhang$driver == 1)))
num_genes = c(length(bone_df$gene), length(zhang$gene))
prop.test(num_drivers, num_genes, alternative = "greater")

#what if we use a different definition of driver? Not zhang definition but hyperBonf
num_drivers = c(length(which(bone_df$hyper_bonf<0.05)), length(which(zhang$hyper_bonf <0.05)))
num_genes = c(length(bone_df$gene), length(zhang$gene))
prop.test(num_drivers, num_genes, alternative = "greater")
#same but for bone genes and NOT bone genes
num_drivers = c(length(which(bone_df$hyper_bonf<0.05)), length(which(not_bone_df$hyper_bonf <0.05)))
num_genes = c(length(bone_df$gene), length(not_bone_df$gene))
prop.test(num_drivers, num_genes, alternative = "greater")
##

#bone gene AND drivers / bone genes
#not bone genes AND drivers / not bone genes
m <- length(bone_df$gene)       # Genes IN bone
n <- length(not_bone_df$gene)       # Genes NOT IN bone
k <- length(which(zhang$driver==1))       # Gene hits, that is, drivers
x <- length(which(bone_df$driver==1))  # Genes both IN bone and drivers
phyper(x-1, m, n, k, log = FALSE,lower.tail = FALSE)##x or x-1?
#
m <- length(bone_df$gene)       # Genes IN bone
n <- length(not_bone_df$gene)       # Genes NOT IN bone
k <- length(which(zhang$hyper_bonf<0.05))       # Gene hits, that is, drivers
x <- length(which(bone_df$hyper_bonf<0.05))  # Genes both IN bone and drivers
phyper(x, m, n, k, log = FALSE,lower.tail = FALSE)##x or x-1?
##
#are key drivers more likely to be highly connected?
key_driver = zhang[which(zhang$hyper_bonf<0.05),]
not_key_driver = zhang[which(zhang$hyper_bonf>=0.05),]

wilcox.test(key_driver$degree, not_key_driver$degree, alternative = "greater")
wilcox.test(key_driver$num_neib, not_key_driver$num_neib, alternative = "greater")
##

#look at module membership vs key driver
#WHY?
#module membership for each genes final module
mainModMem = geneModMemAnnot
colnames(mainModMem)[9:44] = unlist(strsplit(colnames(geneModMemAnnot[,c(9:44)]),split = "ME"))[seq(2,72,by = 2)]
#remove grey genes
mainModMem = mainModMem[-which(mainModMem$color=="grey"),]
mainModMem$MM = NA

for(i in 1:nrow(mainModMem)){
  cx = which(colnames(mainModMem)==mainModMem$color[i])
  mainModMem$MM[i] = mainModMem[i,cx]
}
mainModMem = merge(mainModMem,zhang, by="gene")
#do drivers have greater modMem?
driver = mainModMem[which(mainModMem$driver==1),]
not_driver = mainModMem[which(mainModMem$driver!=1),]

wilcox.test(abs(driver$MM), abs(not_driver$MM), alternative = "less")#actually seems link drivers have lower MM??

driver = mainModMem[which(mainModMem$hyper_bonf<0.05),]
not_driver = mainModMem[which(mainModMem$hyper_bonf>=0.05),]

wilcox.test(abs(driver$MM), abs(not_driver$MM), alternative = "greater")#drivers here have higher MM??

##
##
##
##
#are drivers (zhang) more likely to be colocalizing genes
driver = zhang[which(zhang$driver==1),]
not_driver = zhang[which(zhang$driver!=1),]

m <- length(driver$gene)       # driver genes
n <- length(not_driver$gene)       # non driver
k <- length(which(zhang$coloc_H4>=0.75))       # Gene hits, that is, coloc
x <- length(which(driver$coloc_H4>=0.75))  # Genes both IN driver and coloc
phyper(x-1, m, n, k, log = FALSE,lower.tail = FALSE)##x or x-1?

#are coloc genes more likely to be drivers?
coloc = zhang[which(zhang$coloc_H4>=0.75),]
#not_coloc = zhang[which(zhang$coloc_H4<0.75),]

m <- length(coloc$gene)       # coloc genes
n <- length(zhang$gene) - m       # non coloc
k <- length(which(zhang$driver == 1))       # Gene hits, that is, driver
x <- length(which(coloc$driver==1))  # Genes both IN driver and coloc
phyper(x-1, m, n, k, log = FALSE,lower.tail = FALSE)##x or x-1?

##
x = zhang
#x = x[-which(x$num_neib<2),]

bone_df = x[which(tolower(x$gene) %in% tolower(superset)),]
bone_df$bone = "bone gene"
bone_df$driverx = NA
bone_df[which(bone_df$driver==1),"driverx"] = "driver"
bone_df[which(bone_df$driver==0),"driverx"] = "not_driver"

bone_df$driverz = NA
bone_df[which(bone_df$hyper<0.05),"driverz"] = "driver"
bone_df[which(is.na(bone_df$driverz)),"driverz"] = "not_driver"

not_bone_df = x[which(tolower(x$gene) %in% tolower(superset) == FALSE),]
not_bone_df$bone = "not_bone gene"

not_bone_df$driverx = NA
not_bone_df[which(not_bone_df$driver==1),"driverx"] = "driver"
not_bone_df[which(not_bone_df$driver==0),"driverx"] = "not_driver"

not_bone_df$driverz = NA
not_bone_df[which(not_bone_df$hyper<0.05),"driverz"] = "driver"
not_bone_df[which(is.na(not_bone_df$driverz)),"driverz"] = "not_driver"



df = rbind(bone_df, not_bone_df)
p = ggplot(df,aes(x=num_neib, fill=bone)) + geom_density(position="identity", alpha=0.4, mapping = (aes(y = ..scaled..))) +scale_x_continuous(name = "Number of neighbors",
                                                                                                         breaks = seq(0, 31, 5),
                                                                                                         limits=c(0, 31)) +
  scale_y_continuous(name = "Density") +
  ggtitle("Density plot of number of neighbors") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, family = "Tahoma", face = "bold", hjust=0.5),
        text = element_text(size = 14, family = "Tahoma")) +
  scale_fill_brewer(palette="Dark2", labels=c("Bone genes","Not bone genes"),name="Bone gene status")

compare_means(num_neib ~ bone,  data = df, method = "wilcox.test")


p + annotate("text",x=25, y=.95, label="P=7.1e-09", size=6) 
##stack vs identity?

###stacked bar plot, percentage of drivers bone and not bone
p <- ggplot(data = df,
            mapping = aes(x = bone, fill = driverx))
p + geom_bar(position = position_fill(reverse=T), alpha = 0.8) + ggtitle("Fraction of drivers for bone and non-bone genes") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma")) +
  scale_fill_brewer(palette="Dark2") + labs(fill="Driver Status")



#"p=6.4e-06"
