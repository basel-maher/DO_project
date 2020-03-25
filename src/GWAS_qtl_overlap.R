#BMD GWAS overlap
set.seed(8675309)
library(GenomicRanges)
library(Homo.sapiens)

###
#get mouse QTL loci and identify syntenic regions.
#get qtl
qtl_loc = read.csv("./results/flat/qtl_loc", stringsAsFactors = FALSE)

#convert to BED3 format while getting max interval size in each CI
bed = as.data.frame(matrix(ncol=3, nrow=length(unique(qtl_loc$locus))))
colnames(bed) = c("chr","start","end")

for(i in 1:length(unique(qtl_loc$locus))){
  sub = subset(qtl_loc, locus == unique(qtl_loc$locus)[i])
  start = min(sub$ci_lo)*1000000
  end = max(sub$ci_hi)*1000000
  chrom = paste0("chr",unique(sub$chr))
  bed[i,] = c(chrom, start, end)
}


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



colnames(lifted) = c("hg19_chr","hg19_start","hg19_end")


lifted$locus = 1:11

#merge with qtl loci
lifted = merge(lifted, qtl_loc, by="locus")

#granges format
lifted_mouse_loci = paste0(lifted$hg19_chr,":",lifted$hg19_start, "-",lifted$hg19_end) 
#GRanges
lifted_mouse_loci = as(lifted_mouse_loci, "GRanges")



##load in the Morris associations
#Morris GWAS variants that exceen genome-wide significance threshold (P.NI)

#get ALL Morris GWAS variants
#downloaded from GEFOS (March 3 2020)
#http://www.gefos.org/?q=content/data-release-2018
morris_lead_snps = read.csv("./data/Morrisetal2018.NatGen.SumStats/Morris_eBMD_conditionally_ind_snps.csv", header=T, stringsAsFactors = F)
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

########calculate how many of the 11 loci's syntenic regions would overlap just by chance, 
########by picking 1000 human genomic intervals of the same size distribution as for the 11 human syntenic regions
chrom_lengths = seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)[1:23]
chr = names(chrom_lengths)

lifted_mouse_loci = unique(lifted_mouse_loci)

interval_widths = lifted_mouse_loci@ranges@width


##
sim_frame = as.data.frame(matrix(nrow=1000, ncol=11))
for(i in 1:ncol(sim_frame)){
  chr_sample = sample(x=chr, size = 1000, replace=T, prob = chrom_lengths)
  for(j in 1:1000){
    x = chrom_lengths[chr_sample[j]] - interval_widths[i]
    interval_start = sample(x=1:x,size = 1)

    interval_end = interval_start + interval_widths[i] -1
    
    
    sim_frame[j,i] = paste0(chr_sample[j], ":", interval_start, "-", interval_end)
    print(paste0(i,":",j))
  }
  
}

sim_o = c()
for(i in 1:nrow(sim_frame)){
  g1 = as(sim_frame[i,1], "GRanges")
  g2 = as(sim_frame[i,2], "GRanges")
  g3 = as(sim_frame[i,3], "GRanges")
  g4 = as(sim_frame[i,4], "GRanges")
  g5 = as(sim_frame[i,5], "GRanges")
  g6 = as(sim_frame[i,6], "GRanges") 
  g7 = as(sim_frame[i,7], "GRanges")
  g8 = as(sim_frame[i,8], "GRanges")
  g9 = as(sim_frame[i,9], "GRanges")
  g10 = as(sim_frame[i,10], "GRanges")
  g11 = as(sim_frame[i,11], "GRanges")
  
  g = unlist(GRangesList(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11))
  
  sim_overlaps = GenomicRanges::findOverlaps(query = pos_human_grange, subject = g)
  sim_o = append(sim_o, length(unique(sim_overlaps@to)))
  print(i)
}
quantile(sim_o,probs=0.95)
