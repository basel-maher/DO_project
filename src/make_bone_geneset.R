library(tidyverse)
library(biomaRt)
###

#convert human to mouse
#use mgi homolog list
#http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt (Mar-14-2020)

homology = read.table("./data/mgi_homologs.txt",sep = "\t", header = T)

convertHumantoMouse = function(x){
  mouse=c()
  for(i in 1:length(x)){
    id = homology[which((homology$Common.Organism.Name == "human") & tolower(homology$Symbol) == tolower(x[i])),"HomoloGene.ID"]
    mus = homology[which((homology$Common.Organism.Name == "mouse, laboratory") & homology$HomoloGene.ID == id),"Symbol"]
    mouse=append(mouse,mus)
  }
  mouse=unique(mouse)
  return(mouse)
  
}


###


#Using AmiGO2, downloaded GO terms for the following terms (osteo, bone, ossif). Accessed 7/28/19
#Used filters: is_obsolete:False and idspace: GO
bone_terms = read.delim("./data/GO_term_bone.txt", stringsAsFactors = FALSE, header = FALSE)
#trim to exclude some terms

ex = c("monocyte","megakaryocyte","hair","kidney","neuro","ureter","B cell","tolerance","tendon","muscle","heart","cardio","beak","nephric","tooth","chemotaxis","hemopoiesis","amniotic","wishful")
bone_terms = bone_terms[-(grep(pattern = paste(ex,collapse = "|"), x=bone_terms$V3,ignore.case = TRUE)),]

bone_terms = bone_terms$V1

#
osteo_terms = read.delim("./data/GO_term_osteo*.txt", stringsAsFactors = FALSE, header = FALSE)
osteo_terms = osteo_terms$V1
#
ossif_terms = read.delim("./data/GO_term_ossif*.txt", stringsAsFactors = FALSE, header = FALSE)
ossif_terms = ossif_terms$V1
#
terms = c(ossif_terms, bone_terms, osteo_terms)
terms = unique(terms)


#gets gene symbol, transcript_id and go_id for all genes annotated with terms, from Mus Ensembl
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations


genes_mus = as.data.frame(matrix(ncol=3))
no_annot=c()
for(i in terms){
  print(i)
  gene.data <- unname(getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'go_id'),
                    filters = 'go', values = i, mart = ensembl))
  
  if(nrow(gene.data) >0){
    colnames(gene.data) = colnames(genes_mus)
    genes_mus = rbind(genes_mus,gene.data)
  }else{
    no_annot = append(no_annot, i)
  }
  
}

#same for human
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

genes_hum = as.data.frame(matrix(ncol=3))
no_annot=c()
for(i in terms){
  print(i)
  gene.data <- unname(getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'go_id'),
                            filters = 'go', values = i, mart = ensembl))
  
  if(nrow(gene.data) >0){
    colnames(gene.data) = colnames(genes_hum)
    genes_hum = rbind(genes_hum,gene.data)
  }else{
    no_annot = append(no_annot, i)
  }
  
}
genes_hum = genes_hum[-1,1]
genes_hum = unique(genes_hum)
mussed_human_genes = convertHumantoMouse(genes_hum)

genes_mus= genes_mus[-1,1]
genes_mus = unique(genes_mus)

genes = c(mussed_human_genes, genes_mus)

genes=tolower(genes)
unq_genes = unique(genes)

 
#add MGI genes. manually downloaded osteoporosis, bone mineral density, osteoblast clast and cyte. human and mouse genes
mgi = read.delim("./data/MGIhdpQuery_markers_20190728_224719.txt",stringsAsFactors = FALSE)

mgi_mouse = mgi[which(mgi$Organism=="mouse"),]
mgi_mouse = mgi_mouse$Gene.Symbol
#remove genes with "("
mgi_mouse = mgi_mouse[-grep("\\(", mgi_mouse)]
mgi_mouse = mgi_mouse[-which(mgi_mouse == "917M")]


mgi_hum = mgi[which(mgi$Organism=="human"),]
mgi_hum = mgi_hum$Gene.Symbol
mussed_human_genes = convertHumantoMouse(mgi_hum)

mgi = c(mgi_mouse, mussed_human_genes)

mgi = tolower(mgi)
mgi = unique(mgi)

#superduperset = append(superduperset,mgi)
superduperset = append(mgi, unq_genes)



superduperset = na.omit(superduperset)

superduperset = unique(superduperset)

superduperset = superduperset[-grep("\\.", superduperset)]

write.table(superduperset, "./results/flat/superduperset_sansGWAS.txt", sep = "\t", col.names = FALSE, row.names=FALSE, quote=FALSE)

###

###########Add IMPC genes that have a nominally significant (0.05), weight corrected genotype effect
# impc = read.delim("./data/IMPC_BMD_Results.csv",stringsAsFactors = FALSE, sep=",")
# impc = impc[which(impc$genotype_p_value <0.05),]
# impc = impc$Gene
# impc = unique(impc)
# 
# superduperset = append(superduperset, impc)
# 
# superduperset = na.omit(superduperset)
# 
# superduperset = unique(superduperset)
# 
# superduperset = superduperset[-grep("\\.", superduperset)]
# 
# ##
# write.table(superduperset,"./results/flat/superduperset_GO_MGI_IMPC.txt", sep = "\t", col.names = FALSE, row.names=FALSE, quote=FALSE)





##using the RNA-seq counts (will be used in WGCNA scripts), get the "bone genes" in the superset that are expressed in our samples AFTER filtration. So basically only use the "known bone genes" 
#that were used in the networks analysis

superset = read.delim("./results/flat/superduperset_sansGWAS.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]

# read in the RNA-seq processed counts file
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)


#
annot_file = read.csv("~/Documents/projects/DO_project/results/flat/annot_file.csv")
annot_file = annot_file[,c(1,2)]

counts = counts[which(rownames(counts) %in% annot_file$Gene.ID),]
#

#find and remove features that have fewer than 10 reads in more than 90% (173) of samples 
x=c()
for(i in 1:nrow(counts)){
  if(sum(counts[i,]<10)>=173){
    print(i)
    x = append(x,i)
  }
}

#311 genes removed
counts = counts[-x,]


x = rownames(counts)
x = annot_file[which(annot_file$Gene.ID %in% x),"Gene.Name"]

superset_in_networks = unique(superset[which(tolower(superset) %in% tolower(x))])

write.table(superset_in_networks, "./results/flat/superset_in_networks.txt", sep = "\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
