#BMD overlap
#get all genes around associations, look at overlap with human BMD associated genes.


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





#human BMD genes

############################
#make set from gwascatalog (All associations V1.0.2, accessed 8/22/2019 )
g_catalog = read_tsv("./data/gwas_catalog_v1.0.2-associations_e96_r2019-07-30.tsv")#
terms = c("bone mineral density")
x = g_catalog[grep(pattern = paste(terms,collapse = "|"), g_catalog$`DISEASE/TRAIT`,ignore.case = TRUE),]#filter by terms
# 
# #remove terms that include lead, medication, graft and arthritis, chemotherapy, asthma and alcohol
terms = c("lead","medication", "graft","arthritis","chemotherapy","asthma", "alcohol")
x = x[-(grep(pattern = paste(terms,collapse = "|"), x$`DISEASE/TRAIT`,ignore.case = TRUE)),]
# 
# #
catalog_genes = x$`REPORTED GENE(S)` #get genes
catalog_genes = na.omit(catalog_genes) #remove NAs
#catalog_genes = catalog_genes[-which(catalog_genes=="intergenic")] #remove the term intergenic
catalog_genes = unlist(str_split(catalog_genes,",")) #split elements with multiple genes separated by a comma
catalog_genes = trimws(catalog_genes,"both") #trim whitespace
# 
catalog_genes = tolower(catalog_genes)
catalog_genes = unique(catalog_genes)
# 
mussed_human_genes = convertHumanGeneList(catalog_genes)
mussed_human_genes = unique(mussed_human_genes)

#get qtl list, passed threshold
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

#remove FFP and soleus weight
qtl_norm = qtl_norm[-c(1:3),]

#get genes in CI of qtl
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")

qtl_genes = c()

for(i in 1:nrow(qtl_norm)){
  qtl_genes = c()
  
  #print(i)
  chr = qtl_norm$chr[i]
  start = qtl_norm$ci_lo[i]
  end = qtl_norm$ci_hi[i]
  
  genes_locus <- query_genes(chr, start, end)
  
  qtl_genes = append(qtl_genes, genes_locus$Name)
  if(any(qtl_genes %in% mussed_human_genes)){
    print(qtl_norm$lodcolumn[i])
    print(qtl_genes[which(qtl_genes %in% mussed_human_genes)]
)
  }
}

qtl_genes = unique(qtl_genes)

qtl_genes[which(qtl_genes %in% mussed_human_genes)]





#do same but only for BMD mouse associations

qtl_genes_bmd <- query_genes(8, 102.7427, 104.3691)


qtl_genes_bmd = unique(qtl_genes_bmd$Name)

qtl_genes_bmd[which(qtl_genes_bmd %in% mussed_human_genes)]



# 