library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
###############################################################################################
# integrating human and mouse gwas results
# read in human gwas snps
snps<-read.table("./data/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",header=T)

# read in human snps found in region syntenic with chr. 1 association in the mouse
h = read.csv("./results/flat/gwas_qtl_overlaps.csv")

h1<- h[which(h$locus == 1),]
par(mfrow=c(1,3))

plot(h1$BP,-log10(h1$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr1", ranges = IRanges(start = min(h1$BP), end= max(h1$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}



#2

h2<- h[which(h$locus == 2),]

plot(h2$BP,-log10(h2$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr15", ranges = IRanges(start = min(h2$BP), end= max(h2$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}




#loc 3

h3<- h[which(h$locus == 3),]

plot(h3$BP,-log10(h3$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr20", ranges = IRanges(start = min(h3$BP), end= max(h3$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}


#loc 4

h4<- h[which(h$locus == 4),]

plot(h4$BP,-log10(h4$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr3", ranges = IRanges(start = min(h4$BP), end= max(h4$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}


####loc5


h5<- h[which(h$locus == 5),]

plot(h5$BP,-log10(h5$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr1", ranges = IRanges(start = min(h5$BP), end= max(h5$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}




####loc6


h6<- h[which(h$locus == 6),]

plot(h6$BP,-log10(h6$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr1", ranges = IRanges(start = min(h6$BP), end= max(h6$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}




####loc7


h7<- h[which(h$locus == 7),]

plot(h7$BP,-log10(h7$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,9), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr16", ranges = IRanges(start = min(h7$BP), end= max(h7$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}




####loc8


h8<- h[which(h$locus == 8),]

plot(h8$BP,-log10(h8$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-5,185), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr6", ranges = IRanges(start = min(h8$BP), end= max(h8$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}

####loc8


h8<- h[which(h$locus == 8),]

plot(h8$BP,-log10(h8$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-5,185), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr6", ranges = IRanges(start = min(h8$BP), end= max(h8$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}




####loc9


h9<- h[which(h$locus == 9),]

plot(h9$BP,-log10(h9$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-5,185), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr3", ranges = IRanges(start = min(h9$BP), end= max(h9$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}





####loc10


h10<- h[which(h$locus == 10),]

plot(h10$BP,-log10(h10$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-5,185), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chrX", ranges = IRanges(start = min(h10$BP), end= max(h10$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}


####loc11


h11<- h[which(h$locus == 11),]

plot(h11$BP,-log10(h11$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-5,185), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chrX", ranges = IRanges(start = min(h11$BP), end= max(h11$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=1, 
       col='red',srt=90)
  y0 = y0-0.01
  y1=y1-0.01
}

