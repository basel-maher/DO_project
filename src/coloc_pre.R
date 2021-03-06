#prepare BANs for colocalization and analyze coloc output

options(stringsAsFactors = FALSE)
library(biomaRt)
library(org.Hs.eg.db)


###############################################FUCNTION TO CONVERT MOUSE GENES TO HUMAN HOMOLOGS#####################################
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
###############################################

#Need to prepare genes for colocalization.

#First, get BANs
all = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4_REV.csv", stringsAsFactors = F)
all_f = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4_REV.csv", stringsAsFactors = F)
all_m = read.csv("./results/flat/key_driver_analysis_MALES_sft5_REV.csv", stringsAsFactors = F)


#Analysis
BANs = c(all[which(all$hyper<=0.05), "gene"], all_m[which(all_m$hyper<=0.05), "gene"], all_f[which(all_f$hyper<=0.05), "gene"])
BANs = unique(BANs)

#1465 BANs
#1,370 after comma fix


#get BAN human homologs
BAN_human = convertMousetoHuman(BANs)
#1256
#1,179 after comma fix
#get homologs genomic position
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org",path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")

#need to use entrez gene ids because with ensembl, some gene locations wont show up, as the gene names we have are aliases
gene_pos = getBM(attributes = c("hgnc_symbol","chromosome_name", "start_position","end_position", "ensembl_gene_id", "gene_biotype"),
                 filters = "entrezgene_id",
                 values = BAN_human,
                 mart = mart)


gene_pos = gene_pos[-which(gene_pos$chromosome_name %in% c(1:22,"X")==FALSE),]
gene_pos = gene_pos[-which(gene_pos$hgnc_symbol == ""),] #remove genes with no name 
BAN_human_grange = paste0("chr",gene_pos$chromosome_name, ":", gene_pos$start_position, "-", gene_pos$end_position)
BAN_human_grange = as(BAN_human_grange, "GRanges")
#1251 BANs have human homologs
#1,173 after comma removal
###################################################################
###################################################################
###################################################################
###################################################################
#get eBMD GWAS loci
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
overlaps = GenomicRanges::findOverlaps(query = morris_pos, subject = BAN_human_grange)



#734 BANs overlap morris eBMD loci
#684 after comma removal
homologs_win_1_morris = unique(gene_pos$hgnc_symbol[overlaps@to])
homologs_win_1_ens_morris = unique(gene_pos$ensembl_gene_id[overlaps@to])

write.table(homologs_win_1_morris, file = "./results/flat/homologs_within1mbp_morris_REV2.txt",quote = F,row.names = F,col.names = F)

#format lead snps for colocalization  
snp_gene_df = morris_lead_snps[overlaps@from,]
snp_gene_df$overlap_gene = gene_pos$hgnc_symbol[overlaps@to]
snp_gene_df$overlap_ens = gene_pos$ensembl_gene_id[overlaps@to]
snp_gene_df$overlap_start = gene_pos$start_position[overlaps@to]
snp_gene_df$overlap_end = gene_pos$end_position[overlaps@to]

snp_gene_df = snp_gene_df[,c("SNPID", "RSID","CHR","BP","EA","NEA","MAF","P.NI","N","overlap_gene","overlap_ens")]


snp_gene_df$snpid_19 = paste0(gsub(x=snp_gene_df$SNPID,pattern = ":", replacement = "_" ),"_","b37")


#contains BANs within 1 mbp of morris lead snps, and the snp they are in proximity of
write.table(snp_gene_df[,c("SNPID","snpid_19", "RSID","CHR","BP","overlap_ens","overlap_gene")], file = "./results/flat/coloc/morris_lead_BAN_overlaps_REV2.txt",quote = F,row.names = F,col.names = F)

##############
##############
##############
#Repeat with Estrada lead SNPs. Got lead SNPs from publication, and used GTEx v7 lookup file to convert RSIDs to position for colocalization.

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
overlaps = GenomicRanges::findOverlaps(query = estrada_pos, subject = BAN_human_grange)

length(unique(overlaps@to))
length(unique(overlaps@from))

#105 BANs overlap estrada bmd loci
#99 after comma fix
homologs_win_1_estrada = unique(gene_pos$hgnc_symbol[overlaps@to])
homologs_win_1_ens_estrada = unique(gene_pos$ensembl_gene_id[overlaps@to])

write.table(homologs_win_1_estrada, file = "./results/flat/homologs_within1_estrada_REV2.txt",quote = F,row.names = F,col.names = F)

#format lead snps for colocalization  
snp_gene_df = estrada_lead_snps[overlaps@from,]
snp_gene_df$overlap_gene = gene_pos$hgnc_symbol[overlaps@to]
snp_gene_df$overlap_ens = gene_pos$ensembl_gene_id[overlaps@to]
snp_gene_df$overlap_start = gene_pos$start_position[overlaps@to]
snp_gene_df$overlap_end = gene_pos$end_position[overlaps@to]

#get chrom number, use that to create id_19
snp_gene_df$chr2 = sapply(strsplit(snp_gene_df$chrom,"chr"),"[",2)

snp_gene_df$snpid_19 = paste0(snp_gene_df$chr2,"_",snp_gene_df$BP,"_",snp_gene_df$A1,"_",snp_gene_df$A2,"_","b37")


#contains BANs within 1 mbp of estrada lead snps, and the snp they are in proximity of
write.table(snp_gene_df[,c("RSID","chr2","BP","snpid_19","overlap_ens","overlap_gene")], file = "./results/flat/coloc/estrada_lead_BAN_overlaps_REV2.txt",quote = F,row.names = F,col.names = F)

homologs = unique(c(homologs_win_1_estrada, homologs_win_1_morris))
length(homologs)
#738 BANs within 1mbp of GWAS loci
#688 after comma fix
##convert FNBMD and LSBMD summary stats to RSID
#DO THIS ON CLUSTER
#  lsbmd = read.table("data/GEFOS/GEFOS2_LSBMD_POOLED_GC.txt", header = T, stringsAsFactors = F)
#  fnbmd = read.table("data/GEFOS/GEFOS2_FNBMD_POOLED_GC.txt", header = T, stringsAsFactors = F)
#  
# lookup =  read.table("data/GEFOS/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.lookup_table.txt", header = T, stringsAsFactors = F)
# 
# lsbmd_merged = merge(lsbmd, lookup, by.x = "MarkerName", by.y = "rs_id_dbSNP147_GRCh37p13", all.x=T)
# fnbmd_merged = merge(fnbmd, lookup, by.x = "MarkerName", by.y = "rs_id_dbSNP147_GRCh37p13", all.x=T)
# write.table(lsbmd_merged, file = "lsbmd_merged", row.names = F, quote=F, sep="\t")
# write.table(fnbmd_merged, file = "fnbmd_merged", row.names = F, quote=F, sep="\t")