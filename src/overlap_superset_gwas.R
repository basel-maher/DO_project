#how many known bone genes overlap with GWAS loci, and are known bone genes more likely to overlap with gwas loci than a random set of genes?

#how many known bone genes overlap with GWAS loci

#1-convert superset from mouse to human
superset = read.delim("./results/flat/superset_in_networks.txt", stringsAsFactors = FALSE, header = FALSE)

require(data.table)
MGI_markers = data.table::fread("data/MGI_mouse_genetic_markers_11920.rpt",sep= "\t" , header=T)
MGI_markers = MGI_markers[-which(MGI_markers$Chr=="UN"),]
MGI_markers$`Marker Synonyms (pipe-separated)` = gsub("\\|", " ", MGI_markers$`Marker Synonyms (pipe-separated)`)
MGI_markers$name = tolower(MGI_markers$`Marker Symbol`)

superset = merge(superset, MGI_markers, by.x="V1", by.y="name", all.x=T)

#which(duplicated(superset$V1)) ##no dups anymore
#superset = superset[-1356,]

zz = which(is.na(superset$`MGI Accession ID`) == FALSE)

superset$MGI_NAME[zz] = superset$`Marker Symbol`[zz]

zz = which(is.na(superset$`MGI Accession ID`))
require(Hmisc)
superset$MGI_NAME[zz] = capitalize(superset$V1[zz])

#superset[which(is.na(superset$`MGI Accession ID`)),"MGI Accession ID"] = c("MGI:2140364","MGI:3575190", "MGI:3575247", "MGI:3575248", "MGI:5573174","MGI:1915720","MGI:3819962", "MGI:5633762","MGI:3812132")

#superset[which(superset$MGI_NAME == "Adprhl2"),"Feature Type"] = "protein coding gene"
superset[which(superset$MGI_NAME == "Adprhl2"),] = c("adprhl2","MGI:2140364","4","60.36","126316351","126321703","-","Adprs","O", "ADP-ribosylserine hydrolase","Gene","protein coding gene","Adprhl2 Arh3","Adprhl2")
#superset[which(superset$MGI_NAME == "Impad1"),"Feature Type"] = "protein coding gene"
superset[which(superset$MGI_NAME == "Impad1"),] = c("impad1","MGI:1915720","4","2.56","4762484","4793306","-","Bpnt2","O","3'(2'), 5'-bisphosphate nucleotidase 2","Gene","protein coding gene","gPAPP 1110001C20Rik Jaws Impad1","Impad1")
#superset[which(superset$MGI_NAME %in% c("Hlb324b",  "Hlb328",   "Hlb330",   "Hydro",   "M1665asr", "Rdns","Shsn")),"Feature Type"] = "heritable phenotypic marker"

#superset_pruned = superset[which(superset$`Feature Type` %in% c("complex/cluster/region","heritable phenotypic marker","QTL","unclassified other genome feature")==FALSE),14]
superset=superset[,c(14,2)]

##
homology = read.table("./data/mgi_homologs.txt",sep = "\t", header = T,stringsAsFactors = FALSE)#MGI homology table

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


human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
mouse = useMart("ensembl",dataset="mmusculus_gene_ensembl")

genes = getLDS(attributes = c("mgi_id"),filters="mgi_id", values=superset$`MGI Accession ID`, mart=mouse, attributesL =c("entrezgene_id"),martL=human, uniqueRows=T)
x = superset[which(superset$`MGI Accession ID` %in% genes$MGI.ID == FALSE),]

x$entrez = NA
for(i in 1:nrow(x)){
  id = homology[which((homology$Common.Organism.Name == "mouse, laboratory") & homology$Mouse.MGI.ID == x[i,2]),"HomoloGene.ID"]
  hum = homology[which((homology$Common.Organism.Name == "human") & homology$HomoloGene.ID == id),"EntrezGene.ID"]
  
  if(length(hum) >0){
    x$entrez[i] = hum
  }
}


x[which(is.na(x$entrez)),]

x = x[-which(is.na(x$entrez)),]

colnames(genes) = colnames(x)[c(2,3)]
genes = rbind(genes, x[,c(2,3)])
genes_entrez = unique(genes$entrez)

x = genes[-which(duplicated(genes$`MGI Accession ID`)),]
#fix. too many entrez IDs and not all genes in superset are carried forward
genes_entrez = (unique(x$entrez))
########
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

gene_pos = getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position","end_position", "ensembl_gene_id", "gene_biotype"),
                 filters = "entrezgene_id",
                 values = genes_entrez,
                 mart = mart)


gene_pos = gene_pos[-which(gene_pos$chromosome_name %in% c(1:22,"X")==FALSE),]
gene_pos = gene_pos[-which(gene_pos$hgnc_symbol == ""),] #remove genes with no name 
gene_human_grange = paste0("chr",gene_pos$chromosome_name, ":", gene_pos$start_position, "-", gene_pos$end_position)
gene_human_grange = as(gene_human_grange, "GRanges")



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
morris_pos = paste0(bp_human$chr,":",bp_human$start,"-",bp_human$end)
morris_pos = unique(morris_pos)
morris_pos = as(morris_pos, "GRanges")


#overlap BANs positions with morris lead snps +/- 1Mbp to find BANs that are within 1 Mbp +/- morris gwas lead snps
overlaps = GenomicRanges::findOverlaps(query = morris_pos, subject = gene_human_grange)



#769 known bone genes overlap morris eBMD loci
homologs_win_1_morris = unique(gene_pos$hgnc_symbol[overlaps@to])
homologs_win_1_ens_morris = unique(gene_pos$ensembl_gene_id[overlaps@to])

#get FNBMD and LSBMD GWAS loci
estrada_lead_snps = read.table("./data/GEFOS/lead_snps_pos", header=F, stringsAsFactors = F)
colnames(estrada_lead_snps) = c("chrom","start","end","A1","A2","RSID")
estrada_lead_snps$BP = estrada_lead_snps$start

#make into GRanges format
bp_human = estrada_lead_snps
bp_human$start = bp_human$BP - 1000000
bp_human$start[which(bp_human$start <0)] = 1
bp_human$end = bp_human$BP + 1000000
bp_human$chr = estrada_lead_snps$chrom
#convert to GRanges format
pos_human = paste0(bp_human$chr,":",bp_human$start,"-",bp_human$end)

pos = unique(pos_human)

#GRanges
estrada_pos = as(pos, "GRanges")

#overlap
overlaps = GenomicRanges::findOverlaps(query = estrada_pos, subject = gene_human_grange)

#99 overlap
homologs_win_1_ens_estrada = unique(gene_pos$ensembl_gene_id[overlaps@to])

overlaps_both = c(homologs_win_1_ens_morris, homologs_win_1_ens_estrada)
overlaps_both = unique(overlaps_both)


#get all genes in genome
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

x = getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position","end_position", "ensembl_gene_id", "gene_biotype"),mart = mart)
allGenes = x[which(x$gene_biotype == "protein_coding"),]
#allGenes$chromosome_name = paste0("chr",allGenes$chromosome_name)
allGenes = allGenes[-which(allGenes$chromosome_name %in% c(1:22,"X")==FALSE),]
allGenes = unique(allGenes)
allGenes$chromosome_name = paste0("chr",allGenes$chromosome_name)





#how many genes in the genome overlap all bases covered by BANS += 1 mbp?
#combine
combined_pos = c(estrada_pos, morris_pos)
#reduce
reduced_pos = GenomicRanges::reduce(combined_pos)
#total bases covered
sum(width(reduced_pos))
sum(width(reduce(morris_pos)))
sum(width(reduce(estrada_pos)))


allGenes$pos = paste0(allGenes$chromosome_name,":",allGenes$start_position,"-",allGenes$end_position)

allGenes_pos = as(allGenes$pos, "GRanges")
allGenes_pos = unique(allGenes_pos)

overlaps = GenomicRanges::findOverlaps(query = reduced_pos, subject = allGenes_pos)

length(unique(overlaps@to))

ens_prot_overlap = allGenes[overlaps@to,"ensembl_gene_id"]
ens_prot_overlap = unique(ens_prot_overlap)

#10809 protein coding genes overlap gwas


#771 overlap
q=771
k=1270
m=10809
n=20261-10810
phyper(q-1, m,n,k,lower.tail = F)
#2.9e-8



a = 771 # number known genes that are also overlapping with gwas
b = 1270-a # num known genes that are not overlapping with gwas
c = 10122 #which of the 10809 overlapping genes are not known bone genes (length(which(ens_prot_overlap %in% gene_pos$ensembl_gene_id == FALSE)))  #meaning they prot coding genes that overlap gwas and are not known bone genes
d = 8873 #not known bone genes and also not overlapping length(which( (allGenes$ensembl_gene_id %in% ens_prot_overlap == FALSE) & (allGenes$ensembl_gene_id %in% gene_pos$ensembl_gene_id == FALSE) ))

#make contingency table
mat = matrix(nrow=2,ncol=2)
mat[1,] = c(a,c)
mat[2,] = c(b,d)
fisher.test(mat) #two sided (enrichment or depletion)
fisher.test(mat,alternative = "g") #enrichment

x = fisher.test(mat,alternative = "g")
x$p.value

######
a = 738 # BANs that have human homologs and are within 1 mbp of gwas snp
b = 1251-738 #BANs that are human homologs but not within 1 mbp of gwas snp
c = 10809-738 # genes within 1 mbp of gwas snp but arent BANs
d = 20261 - (a+b+c)#not bans and not within 1 mbp

mat = matrix(nrow=2,ncol=2)
mat[1,] = c(a,c)
mat[2,] = c(b,d)
fisher.test(mat) #two sided (enrichment or depletion)
fisher.test(mat,alternative = "g") #enrichment

x = fisher.test(mat,alternative = "g")
x$p.value
#change in rev1.R

