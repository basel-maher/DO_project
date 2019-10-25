library(biomaRt)
library(tidyverse)
###
#from https://www.r-bloggers.com/converting-mouse-to-human-gene-names-with-biomart-package/

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
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
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses mus ensembl annotations

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
mussed_human_genes = convertHumanGeneList(genes_hum)

genes_mus= genes_mus[-1,1]
genes_mus = unique(genes_mus)

genes = c(mussed_human_genes, genes_mus)

genes=tolower(genes)
unq_genes = unique(genes)

#unq_genes = toupper(unq_genes)
############################
#make set from gwascatalog (All associations V1.0.2, accessed 8/22/2019 )
g_catalog = read_tsv("./data/gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv")#
terms = c("bone","bone mineral density", "osteoporo","osteoblast","osteoclast","osteocy")
x = g_catalog[grep(pattern = paste(terms,collapse = "|"), g_catalog$`DISEASE/TRAIT`,ignore.case = TRUE),]#filter by terms

#remove terms that include lead, medication, graft and arthritis, chemotherapy, asthma and alcohol
terms = c("lead","medication", "graft","arthritis","chemotherapy","asthma", "alcohol")
x = x[-(grep(pattern = paste(terms,collapse = "|"), x$`DISEASE/TRAIT`,ignore.case = TRUE)),]

#
catalog_genes = x$`REPORTED GENE(S)` #get genes
catalog_genes = na.omit(catalog_genes) #remove NAs
catalog_genes = catalog_genes[-which(catalog_genes=="intergenic")] #remove the term intergenic
catalog_genes = unlist(str_split(catalog_genes,",")) #split elements with multiple genes separated by a comma
catalog_genes = trimws(catalog_genes,"both") #trim whitespace

catalog_genes = tolower(catalog_genes)
catalog_genes = unique(catalog_genes)

mussed_human_genes = convertHumanGeneList(catalog_genes)
mussed_human_genes = unique(mussed_human_genes)
#
superduperset = append(genes, mussed_human_genes)
superduperset = unique(superduperset)

#add MGI genes. manually downloaded osteoporosis, bone mineral density, osteoblast clast and cyte. human and mouse genes
mgi = read.delim("./data/MGIhdpQuery_markers_20190728_224719.txt",stringsAsFactors = FALSE)

mgi_mouse = mgi[which(mgi$Organism=="mouse"),]
mgi_mouse = mgi_mouse$Gene.Symbol
#remove genes with "("
mgi_mouse = mgi_mouse[-grep("\\(", mgi_mouse)]
mgi_mouse = mgi_mouse[-which(mgi_mouse == "917M")]


mgi_hum = mgi[which(mgi$Organism=="human"),]
mgi_hum = mgi_hum$Gene.Symbol
mussed_human_genes = convertHumanGeneList(mgi_hum)

mgi = c(mgi_mouse, mussed_human_genes)

mgi = tolower(mgi)
mgi = unique(mgi)

superduperset = append(superduperset,mgi)

superduperset = na.omit(superduperset)

superduperset = unique(superduperset)
##
write.table(superduperset,"./results/flat/superduperset.txt", sep = "\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
