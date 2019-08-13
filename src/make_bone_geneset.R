library(biomaRt)

bone_terms = read.delim("./data/GO_term_bone.txt", stringsAsFactors = FALSE, header = FALSE)
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

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl") #uses mus ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with terms

genes = as.data.frame(matrix(ncol=3))
no_annot=c()
for(i in terms){
  print(i)
  gene.data <- unname(getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'go_id'),
                    filters = 'go', values = i, mart = ensembl))
  
  if(nrow(gene.data) >0){
    colnames(gene.data) = colnames(genes)
    genes = rbind(genes,gene.data)
  }else{
    no_annot = append(no_annot, i)
  }
  
}

#getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'go_id'),
#      filters = 'go', values = "GO:0036072", mart = ensembl)

genes = genes[-1,]

unq_genes = unique(genes$V1)

#unq_genes = toupper(unq_genes)
#
#superset from mouse and human gwas
superset = read.csv("~/Desktop/DO_proj/data/networks/human_and_mouse_bone_superset.csv", stringsAsFactors = FALSE, header=FALSE)
superset = superset[,1]
superset = append(superset,"CPE") #for carboxyl peptidase E
superset = unique(superset)
#
superduperset = append(superset, unq_genes)
superduperset = unique(superduperset)

#add MGI genes. manually downloaded osteoporosis, bone mineral density, osteoblast clast and cyte. human and mouse genes
mgi = read.delim("~/Downloads/MGIhdpQuery_markers_20190728_224719.txt",stringsAsFactors = FALSE)
mgi = mgi$Gene.Symbol
#remove genes with "("
mgi = mgi[-grep("\\(", mgi)]
mgi = toupper(mgi)
mgi = unique(mgi)
superduperset = append(superduperset,mgi)
superduperset = superduperset[-c(1444)] #contain \
superduperset = superduperset[-c(1660)] #contain \
superduperset = toupper(superduperset)

superduperset = unique(superduperset)
##
write.table(superduperset,"~/Desktop/superduperset.txt", sep = "\t", col.names = FALSE, row.names=FALSE, quote=FALSE)
