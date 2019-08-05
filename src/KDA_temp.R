library("igraph")
library("bnlearn")
library("parallel")
options(stringsAsFactors = FALSE)

load("./results/Rdata/geneModMemAnnot.RData")
#
bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}

#gene superset
superset = read.delim("~/Desktop/superduperset.txt", stringsAsFactors = FALSE, header = FALSE)
superset = superset[,1]

#BNs learned on high performance computing cluster
bns = list.files("~/Desktop/bn/")
#
x = bn2igraph(red_bn)
subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,2,nodes = "Pparg"))))
plot(subgraph,vertex.label.cex=0.9,edge.width=2, vertex.size=25, margin=-0.1, vertex.label.dist=0.2, vertex.label.degree=-pi)

####

out = list()
counter = 1
for(net in bns){
  print(counter)
  color = strsplit(net,"_")[[1]][2]
  #print(color)
  load(paste0("~/Desktop/bn/",net))
  obj_name = paste0(color,"_bn")
  assign(x = obj_name ,bn)
  
  z = bn2igraph(get(obj_name)) #get "gets" an object from env based on name. otherwise its just a string here
  mod_genes = geneModMemAnnot[which(geneModMemAnnot$color == color),"gene"]
  
  
  kda = as.data.frame(matrix(nrow=length(mod_genes),ncol = 6))
  
  i=0
  
  for(gene in mod_genes){
    
    neighbors = names(unlist(neighborhood(z,nodes = gene,order = 2)))
    num_neib = neighborhood.size(z,nodes = gene,order = 2)
    num_bone_neib = length(which(tolower(neighbors) %in% tolower(superset)))
    in_neib_not_bone = num_neib - num_bone_neib
    not_neib_in_bone = length(which(tolower(mod_genes) %in% tolower(superset) & tolower(mod_genes) %in% tolower(neighbors) == FALSE))
    #not_neib_in_bone =  length(which(tolower(superset[,1]) %in% tolower(neighbors) == FALSE))
    not_neib_not_bone = length(which(tolower(mod_genes) %in% tolower(superset) == FALSE & tolower(mod_genes) %in% tolower(neighbors) == FALSE))
    num_in = length(neighbors(z, gene, mode="in")) #0 "in" arcs signifies a root node
    
    i=i+1
    kda[i,1] = gene
    kda[i,2] = num_neib
    kda[i,3] = num_bone_neib
    kda[i,4] = num_in
    kda[,5] = kda[,3]/kda[,2]
    
    x = matrix(c(num_bone_neib,not_neib_in_bone,in_neib_not_bone,not_neib_not_bone),nrow = 2,ncol = 2)
    xx = fisher.test(x,alternative = "greater")
    kda[i,6] = xx$p.value
    
  }
  
  colnames(kda) = c("gene","num_neib","num_bone_neib","num_in_arcs","ratio_bone_neib","fisher_pval")
  
  kda$bonf_pval = kda$fisher_pval*kda$num_neib
  
  out[[counter]] = kda
  names(out)[counter] = color
  counter = counter+1
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
#annot_file = read.delim("367-FarberDO2_S20.gene_abund.tab",header = TRUE)
#allGenes = unique(annot_file$Gene.Name)
#allGenes = colnames(resid)

allGenes = geneModMemAnnot$gene
perm=0
m=c()
while(perm < 10){
  print(perm)
  #rand_set = sample(allGenes, ncol(blue_genes_exp))
  rand_set = sample(allGenes, length(superset))
  kda_rand = as.data.frame(matrix(nrow=length(mod_genes),ncol = 4))
  i=0
  for(gene in mod_genes){
    neighbors = names(unlist(neighborhood(z,nodes = gene,order = 2)))
    num_neib = neighborhood.size(z,nodes = gene,order = 2)
    num_bone_neib = length(which(tolower(neighbors) %in% tolower(rand_set)))
    num_in = length(neighbors(z, gene, mode="in")) #0 "in" arcs signifies a root node
    i=i+1
    kda_rand[i,1] = gene
    kda_rand[i,2] = num_neib
    kda_rand[i,3] = num_bone_neib
    kda_rand[i,4] = num_in
    
  }
  
  kda_rand[,5] = kda_rand[,3]/kda_rand[,2]
  if(any(kda_rand[,2]==1)&kda_rand[,3]==1))
  #x = kda_rand[-which(kda_rand[,2]==1 & kda_rand[,3]==1),] #FIX
  #m = c(m,as.numeric(x[,5]))
  m = c(m,as.numeric(kda_rand[,5]))
  #m = as.numeric(na.omit(m))
  perm = perm+1
}


#######