library(qtl2)
load("results/Rdata/cross_eqtl.Rdata")
source("./bin/NEO_orig.txt")#https://horvath.genetics.ucla.edu/html/aten/NEO/
#source("bin/CausalityFunctions.txt")
annot_mstrg = read.delim("./data/314-FarberDO2_S6.gene_abund.tab",header = TRUE)
#load("./results/Rdata/trans_eqtl_ALL.Rdata")
mediators = read.csv("~/Desktop/mediators.csv", stringsAsFactors = FALSE)
load("./results/Rdata/cis_eqtl_ALL.Rdata")
#try with mediators where trans-eqtl is within 3 mbp window of mediator cis
mediators_near  = mediators[which(abs(mediators$pos - mediators$trans_orig_pos) <= 1.5),]

#is a transeqtl's expression mediated by its mediators expression? SNP (causal anchor) is original transeqtl (why?)
#makes more sense to use the cis eqtl. info travels cis -> mediator -> target.
#do both
annot_mstrg = annot_mstrg[-which(annot_mstrg$Gene.ID %in% colnames(cross_eqtl$pheno) == FALSE),]

neo_out = as.data.frame(matrix(ncol=37),nrow = nrow(mediators_near))

for(i in 1:nrow(mediators_near)){
  print(i)
  target = mediators_near$target[i]
  mediator = mediators_near$id[i]
  #anchor_chr = mediators_near$trans_orig_chr[i]
  #anchor_pos = mediators_near$trans_orig_pos[i]
  #use actual cis-EQTL markers here
  if(any(which((cis_eqtl_ALL$chr == mediators_near$chr[i]) & (cis_eqtl_ALL$lodcolumn==mediators_near$id[i])))){
    anchor_pos = cis_eqtl_ALL[which((cis_eqtl_ALL$chr == mediators_near$chr[i]) & (cis_eqtl_ALL$lodcolumn==mediators_near$id[i])),"pos"]
  } else {
    next
  }
  anchor_chr = mediators_near$chr[i]
  #anchor_pos = mediators_near$pos[i]
  
  marker = find_marker(cross_eqtl$pmap, chr = anchor_chr,pos = anchor_pos)
  SNP_geno = cross_eqtl$geno[[as.character(anchor_chr)]][,marker[1]]
  
  datSNP = data.frame(SNP_geno)
  
  #recode SNPs
  datSNP[which(datSNP$SNP1==0),] = NA
  datSNP[which(datSNP$SNP1==1),] = 0
  datSNP[which(datSNP$SNP1==2),] = 1
  datSNP[which(datSNP$SNP1==3),] = 2
  
  #make E1 and E2
  id = as.character(annot_mstrg[which(annot_mstrg$Gene.Name == mediator),"Gene.ID"][1])#
  E1 = cross_eqtl$pheno[,id]

  id = as.character(annot_mstrg[which(annot_mstrg$Gene.Name == target),"Gene.ID"][1])#
  E2 = cross_eqtl$pheno[,id]
  
  #combine into 1 df
  datExpr1= data.frame(E1, E2)
  datCombined=data.frame(datExpr1, datSNP)
  
  datCombined = na.omit(datCombined)
  
  #get  and set params
  pm=neo.get.param()
  
  pm$A=1 # E1 corresponds to column 1 in datCombined
  pm$B=2 # E2 corresponds to column 2 in datCombined
  pm$top.N.snps.per.trait=3 # this is the maximum number of anchors for each trait.
  pm$work.hard.at.EO.contingencies.if.not.both.M.B.and.M.A.markers = TRUE
  # this activates the calculation of an edge orienting score, which is explained below.
  pm$do.m1m2.average = TRUE
  
  #single marker analysis
  sma = single.marker.analysis(datCombined,snpcols = "SNP_geno",genecols = "E1",traitcols = "E2",pm = pm)
  
  neo_out[i,c(1:5)] = c(target, mediator, anchor_chr,anchor_pos, marker)
  colnames(neo_out)[1:5] = c("target", "mediator", "anchor_chr","anchor_pos", "marker")
  
  neo_out[i,c(6:37)] = sma
  colnames(neo_out)[6:37] = colnames(sma)
}
#save(neo_out,file="~/Desktop/neo_out.RData")
neo_out = neo_out[-which(is.na(neo_out$target)),]
save(neo_out,file="~/Desktop/neo_out_ciseqtl.RData")

#neo_out_trans = neo_out
#neo_out_cis = neo_out
neo_out_cis_near = neo_out
neo_out_trans_near = neo_out
#



#neo1=NEO.scores(datSNP=datSNP ,A=E1,B=E2, detailed.output=T) 

#####

#get the SNPs around the bglap2 distal eqtl
marker = find_marker(cross_eqtl$pmap, chr = 16,pos = c(24.71717))
marker = find_marker(cross_eqtl$pmap, chr = 3,pos = c(85.52581))

marker = find_marker(cross_eqtl$pmap, chr = 16,interval = c(24.71717-0.25,24.71717+0.25))
marker = find_marker(cross_eqtl$pmap, chr = 10,interval = c(	82.18774-0.25,	82.18774+0.25))#24.71717
marker = find_marker(cross_eqtl$pmap, chr = 10,pos = c(	82.18774))
marker = find_marker(cross_eqtl$pmap, chr = 1,pos = c(100.2614))
SNP1 = cross_eqtl$geno$`3`[,marker[1]]
SNP2 = cross_eqtl$geno$`10`[,marker[2]]
SNP3 = cross_eqtl$geno$`10`[,marker[3]]
SNP4 = cross_eqtl$geno$`10`[,marker[4]]
SNP5 = cross_eqtl$geno$`10`[,marker[5]]
SNP6 = cross_eqtl$geno$`10`[,marker[6]]
SNP7 = cross_eqtl$geno$`10`[,marker[7]]
SNP8 = cross_eqtl$geno$`10`[,marker[8]]
SNP9 = cross_eqtl$geno$`10`[,marker[9]]
SNP10 = cross_eqtl$geno$`10`[,marker[10]]
SNP11 = cross_eqtl$geno$`10`[,marker[11]]
SNP12 = cross_eqtl$geno$`10`[,marker[12]]
SNP13 = cross_eqtl$geno$`10`[,marker[13]]
SNP14 = cross_eqtl$geno$`10`[,marker[14]]
SNP15 = cross_eqtl$geno$`10`[,marker[15]]
SNP16 = cross_eqtl$geno$`10`[,marker[16]]

datSNP = data.frame(SNP1,SNP2,SNP3,SNP4,SNP5,SNP6,SNP7,SNP8,SNP9,SNP10,SNP11,SNP12,SNP13,SNP14,SNP15,SNP16)

datSNP[which(datSNP$SNP1==0),] = NA
datSNP[which(datSNP$SNP1==1),] = 0
datSNP[which(datSNP$SNP1==2),] = 1
datSNP[which(datSNP$SNP1==3),] = 2
#Bglap2

id = as.character(annot_mstrg[which(annot_mstrg$Gene.Name == "Rps3a1"),"Gene.ID"])
E1 = cross_eqtl$pheno[,id]

#
#Trp63
id = as.character(annot_mstrg[which(annot_mstrg$Gene.Name == "Rps3a3"),"Gene.ID"])
E2 = cross_eqtl$pheno[,id]
#
datExpr1= data.frame(E1, E2)
datCombined=data.frame(datExpr1, datSNP)

datCombined = na.omit(datCombined)
pm=neo.get.param()
#
# Here tell NEO that we want to evaluate the edge E1 and E2
pm$A=1 # E1 corresponds to column 1 in datCombined
pm$B=2 # E2 corresponds to column 2 in datCombined
pm$top.N.snps.per.trait=3 # this is the maximum number of anchors for each trait.
pm$work.hard.at.EO.contingencies.if.not.both.M.B.and.M.A.markers = TRUE
# this activates the calculation of an edge orienting score, which is explained below.
pm$do.m1m2.average = TRUE
NEOoutput1 = neo(datCombined,pm=pm)

NEOoutput1$stats[,c("edge", "LEO.NB.OCA", "LEO.NB.CPA", "LEO.NB.MAX", "M1M2.AVG")]
NEOoutput1$stats[NEOoutput1$stats$edge=="E1 -> E2", "M1M2.AVG"]

sma = single.marker.analysis(datCombined,snpcols = "SNP1",genecols = "E1",traitcols = "E2",pm = pm)
sma


