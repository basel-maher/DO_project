#BMD GWAS overlap
library(GenomicRanges)
###
#get mouse QTL loci and identify syntenic regions.
#get qtl
qtl_loc = read.csv("./results/flat/qtl_loc", stringsAsFactors = FALSE)

#convert to BED3 format
bed = as.data.frame(matrix(ncol=3, nrow=nrow(qtl_loc)))
colnames(bed) = c("chr","start","end")                    
bed$chr = qtl_loc$chr
bed$chr = paste0("chr",bed$chr)

bed$start = qtl_loc$ci_lo *1000000
bed$end = qtl_loc$ci_hi*1000000

#convert using UCSC liftOver, from mm10 to hg19
#https://genome.ucsc.edu/cgi-bin/hgLiftOver
#PARAMS USED:

# Minimum ratio of bases that must remap: 	0.10
# 
# BED 4 to BED 6 Options
# Allow multiple output regions: 	off
# Minimum hit size in query: 	100000
# Minimum chain size in target: 	0
# 
# BED 12 Options
# Min ratio of alignment blocks or exons that must map: 	1.00
# If thickStart/thickEnd is not mapped, use the closest mapped base: 	off
#accessed MAR 8 2020

#saved as "./results/flat/lifted_qtl_loci.bed"

##read in the lifted over regions
lifted = read.table("./results/flat/lifted_qtl_loci.bed")

#some loci overlap

colnames(lifted) = c("hg19_chr","hg19_start","hg19_end")

#merge with qtl loci
lifted = cbind(lifted, qtl_loc)

#granges format
lifted_mouse_loci = paste0(lifted$hg19_chr,":",lifted$hg19_start, "-",lifted$hg19_end) 
#GRanges
lifted_mouse_loci = as(lifted_mouse_loci, "GRanges")



##load in the Morris associations
#Morris GWAS variants that exceen genome-wide significance threshold (P.NI)

#get ALL Morris GWAS variants
#downloaded from GEFOS (March 3 2020)
#http://www.gefos.org/?q=content/data-release-2018

morris = read.table("./data/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt", header=T)
#get vars that are genome wide significant
morris = morris[which(morris$P.NI <= 5e-8),]

#create GRanges object
chroms_human = paste0("chr",morris$CHR)
#convert chr23 to chrX
chroms_human[which(chroms_human == "chr23")] = "chrX"
bp_human = morris$BP

#convert to GRanges format
pos_human = paste0(chroms_human,":",bp_human)

#GRanges
pos_human_grange = as(pos_human, "GRanges")


overlaps = GenomicRanges::findOverlaps(query = pos_human_grange, subject = lifted_mouse_loci)

overlaps_human = overlaps@from
overlaps_mouse = overlaps@to  

x = qtl_loc[overlaps_mouse,]

overlaps_merged = cbind(x,morris[overlaps_human,])
write.csv(overlaps_merged, file="./results/flat/gwas_qtl_overlaps.csv",quote = F,row.names = F)

