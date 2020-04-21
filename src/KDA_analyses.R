options(stringsAsFactors = FALSE)
library(topGO)
library(biomaRt)
library(org.Hs.eg.db)
#kda_analyses#
all = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv", stringsAsFactors = F)
all_f = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv", stringsAsFactors = F)
all_m = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv", stringsAsFactors = F)

#bone gene superset
superset = read.delim("./results/flat/superduperset_sansGWAS.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]


#Analysis
genes = c(all[which(all$hyper<=0.05), "gene"], all_m[which(all_m$hyper<=0.05), "gene"], all_f[which(all_f$hyper<=0.05), "gene"])
genes = unique(genes)

length(which(tolower(genes) %in% tolower(superset)))
length(which(tolower(genes) %in% tolower(superset)==FALSE))
##are bone genes more highly connected in the network?
x = all #repeat for males and females

bone_df = x[which(tolower(x$gene) %in% tolower(superset)),]
not_bone_df = x[which(tolower(x$gene) %in% tolower(superset) == FALSE),]


wilcox.test(bone_df$degree, not_bone_df$degree, alternative = "greater")#1.8e-4, 0.0283 males, 1.052e-5 females
wilcox.test(bone_df$num_neib, not_bone_df$num_neib, alternative = "greater")#2.621e-5, 0.005216 males, 5.542e-7 females



###GO analysis
load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")

allGenes = unique(c(combat_annot$Gene.Name, combat_annot_f$Gene.Name, combat_annot_m$Gene.Name))

interesting.genes = not_bone_df$gene

geneList<-factor(as.integer(allGenes %in% interesting.genes)) #If TRUE returns 1 as factor, otherwise 0
names(geneList)<-allGenes
###MF###
GOdata <- new("topGOdata", ontology = "MF", allGenes =geneList,
              annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t1<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t1)
###CC###
GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
              annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t2<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t2)
###BP###
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
result<-getSigGroups(GOdata,test.stat)
t3<-GenTable(GOdata, classic=result, topNodes=length(result@score))
head(t3)
####
t.all = NULL
t.all<-rbind(t1,t2,t3)
t.all$classic<-as.numeric(as.character(t.all$classic))
######
######
######
######
##How many nominal key drivers are located within BMD GWAS loci?

homology = read.table("./data/mgi_homologs.txt",sep = "\t", header = T)

convertMousetoHuman = function(x){
  human=c()
  for(i in 1:length(x)){
    id = homology[which((homology$Common.Organism.Name == "mouse, laboratory") & tolower(homology$Symbol) == tolower(x[i])),"HomoloGene.ID"]
    hum = homology[which((homology$Common.Organism.Name == "human") & homology$HomoloGene.ID == id),"EntrezGene.ID"]
    human=append(human,hum)
  }
  human=unique(human)
  return(human)
  
}


#get nominal genes
#genes object from above

#get homologs
nominal_hum = convertMousetoHuman(genes)

#get homolog location
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

gene_pos = getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position","end_position"),
             filters = "entrezgene_id",
             values = nominal_hum,
             mart = mart)

gene_pos = gene_pos[-which(gene_pos$chromosome_name %in% c(1:22,"X")==FALSE),]

gene_pos_grange = paste0("chr",gene_pos$chromosome_name, ":", gene_pos$start_position, "-", gene_pos$end_position)
gene_pos_grange = as(gene_pos_grange, "GRanges")

#need to use entrez gene ids because with ensembl, some gene locations wont show up, as the gene names we have are aliases

#get BMD GWAS loci
morris_lead_snps = read.csv("./data/Morrisetal2018.NatGen.SumStats/Morris_eBMD_conditionally_ind_snps.csv", header=T, stringsAsFactors = F)

#make into GRanges format
chroms_human = paste0("chr",morris_lead_snps$CHR)
#convert chr23 to chrX
chroms_human[which(chroms_human == "chr23")] = "chrX"
bp_human = as.data.frame(morris_lead_snps$BP)
bp_human$start = bp_human$`morris_lead_snps$BP` - 1000000
bp_human$start[which(bp_human$start <0)] = 1
bp_human$end = bp_human$`morris_lead_snps$BP` + 1000000
bp_human$chr = chroms_human
#convert to GRanges format
pos_human = paste0(bp_human$chr,":",bp_human$start,"-",bp_human$end)

#GRanges
pos_human_grange = as(pos_human, "GRanges")

#overlap
overlaps = GenomicRanges::findOverlaps(query = pos_human_grange, subject = gene_pos_grange)

length(unique(overlaps@to))

#545 genes overlap 680 loci

#permute. how many expected by chance?
allGenes = c(all$gene, all_m$gene, all_f$gene)
allGenes = unique(allGenes)

#because some genes arent in both mouse an human for some reason, get only genes that have a human homolog
allGenesx = allGenes[which(allGenes %in% homology$Symbol)]

allGenesx = homology[which(homology$Symbol %in% allGenesx & homology$Common.Organism.Name == "mouse, laboratory"),"HomoloGene.ID"]

allGenesx = homology[which(homology$HomoloGene.ID %in% allGenesx & homology$Common.Organism.Name == "human"),"HomoloGene.ID"]

allGenesx = homology[which(homology$HomoloGene.ID %in% allGenesx & homology$Common.Organism.Name == "mouse, laboratory"),"Symbol"]



i=0
permx = c()
while(i<1000){
  g = sample(allGenesx,size = 904, replace = F)
  gg  = convertMousetoHuman(g)

  gene_pos = getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position","end_position"),
                   filters = "entrezgene_id",
                   values = gg,
                   mart = mart)
  
  gene_pos = gene_pos[-which(gene_pos$chromosome_name %in% c(1:22,"X")==FALSE),]
  
  gene_pos_grange = paste0("chr",gene_pos$chromosome_name, ":", gene_pos$start_position, "-", gene_pos$end_position)
  gene_pos_grange = as(gene_pos_grange, "GRanges")
  print(length(gene_pos_grange))
  overlaps = GenomicRanges::findOverlaps(query = pos_human_grange, subject = gene_pos_grange)
  
  permx = append(permx,length(unique(overlaps@to)))
  print(i)
  i=i+1
}

quantile(perm,probs=0.95)
save(perm, file="./results/Rdata/perm_kda_in_gwasloci.Rdata")
