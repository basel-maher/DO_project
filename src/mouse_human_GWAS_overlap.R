#BMD overlap
#get all genes around associations, look at overlap with human BMD associated genes.

library(tidyverse)
library(biomaRt)
library(liftOver)
library(GenomicRanges)
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



##################
#Identify syntenic regions with liftOver, between Morris BMD GWAS and our associations

#read in the morris eBMD GWAS data / lead snps from supplementary table 2
#PMCID: PMC6358485
#hg19
morris = read.csv("./data/Morris_eBMD_conditionally_ind_snps.csv", header=T)

#take lead snp positions and convert to GRanges object

#take snps
chroms_human = paste0("chr",morris$CHR)
#convert chr23 to chrX
chroms_human[which(chroms_human == "chr23")] = "chrX"
bp_human = morris$BP

#GRanges format
pos_human = paste0(chroms_human,":",bp_human)

#GRanges
pos_human_grange = as(pos_human, "GRanges")
mcols(pos_human_grange) = pos_human

###liftover from hg19 to mouse mm10
#get chain object
#downloaded from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/ on (Mar 2 2020)
chain <- import.chain("./data/hg19ToMm10.over.chain")

lifted = liftOver(pos_human_grange, chain)
#lifted = unlist(lifted)


#convert our associations, including CI, to GRanges
chr = qtl_loc$chr
chr = paste0("chr",chr)

bp_start = qtl_loc$ci_lo *1000000
bp_end = qtl_loc$ci_hi*1000000
bp_range = paste0(bp_start,"-",bp_end)

pos = paste0(chr,":",bp_range)

mouse = as(pos, "GRanges")

overlaps = findOverlaps(query = lifted, subject = mouse)

overlaps_human = overlaps@from
overlaps_mouse = overlaps@to  

x = qtl_loc[overlaps_mouse,]
y = pos_human_grange@elementMetadata[overlaps_human,]  

overlaps_merged = cbind(x,morris[overlaps_human,])

