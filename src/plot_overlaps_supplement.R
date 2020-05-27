library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
###############################################################################################
# integrating human and mouse gwas results
# read in human gwas snps
#snps<-read.table("./data/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",header=T)

# read in human snps found in region syntenic with chr. 1 association in the mouse
h = read.csv("./results/flat/gwas_qtl_overlaps_allIncNonSig.csv", stringsAsFactors = FALSE)

h1<- h[which(h$locus == 1),]


cairo_pdf(file="~/Desktop/figs/sup1a.pdf", width = 10, height = 7)
plot(h1$BP,-log10(h1$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-3,max(-log10(h1$P.NI))+1), lwd=1)
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

genes = genes[-which(genes$symbol %in% c("LINC01699", "RASAL2-AS1", "MIR3121", "LXH4-AS1", "MIR4424")),]

y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()


#2

h2<- h[which(h$locus == 2),]

cairo_pdf(file="~/Desktop/figs/sup1b.pdf", width = 10, height = 7)
  plot(h2$BP,-log10(h2$P.NI),type="p",
       col= 'red', cex=1.5,
       axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,10), lwd=1)
  abline(h=7.30103)
  
  mapping = as.data.frame(org.Hs.egSYMBOL)
  
  gr = GRanges(seqnames = "chr20", ranges = IRanges(start = min(h2$BP), end= max(h2$BP)))
  
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
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()


#loc 3

h3<- h[which(h$locus == 3),]


cairo_pdf(file="~/Desktop/figs/sup1c.pdf", width = 10, height = 7)
plot(h3$BP,-log10(h3$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-3,max(-log10(h3$P.NI))+1), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr3", ranges = IRanges(start = min(h3$BP), end= max(h3$BP)))

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
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()


#loc 4

h4<- h[which(h$locus == 4),]


cairo_pdf(file="~/Desktop/figs/sup1d.pdf", width = 10, height = 7)
plot(h4$BP,-log10(h4$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-3,max(-log10(h4$P.NI))+1), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr1", ranges = IRanges(start = min(h4$BP), end= max(h4$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

genes = genes[-which(genes$width<10000),]
#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)


y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()



####loc5


h5<- h[which(h$locus == 5),]




cairo_pdf(file="~/Desktop/figs/sup1e.pdf", width = 10, height = 7)
plot(h5$BP,-log10(h5$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-3,8), lwd=1)
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
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()






####loc6


h6<- h[which(h$locus == 6),]



cairo_pdf(file="~/Desktop/figs/sup1f.pdf", width = 10, height = 7)
plot(h6$BP,-log10(h6$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-3,max(-log10(h6$P.NI))+1), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr16", ranges = IRanges(start = min(h6$BP), end= max(h6$BP)))

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
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()




####loc7


h7<- h[which(h$locus == 7),]





cairo_pdf(file="~/Desktop/figs/sup1g.pdf", width = 10, height = 7)
plot(h7$BP,-log10(h7$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-12,max(-log10(h7$P.NI))+1), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr6", ranges = IRanges(start = min(h7$BP), end= max(h7$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)


y0=y1=-8
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()




####loc8


h8<- h[which(h$locus == 8),]





cairo_pdf(file="~/Desktop/figs/sup1h.pdf", width = 10, height = 7)
plot(h8$BP,-log10(h8$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,max(-log10(h8$P.NI))+1), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chr3", ranges = IRanges(start = min(h8$BP), end= max(h8$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap

#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)


y0=y1=-2
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-2,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()



####loc9


h9<- h[which(h$locus == 9),]





cairo_pdf(file="~/Desktop/figs/sup1i.pdf", width = 10, height = 7)
plot(h9$BP,-log10(h9$P.NI),type="p",
     col= 'red', cex=1.5,
     axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,8), lwd=1)
abline(h=7.30103)

mapping = as.data.frame(org.Hs.egSYMBOL)

gr = GRanges(seqnames = "chrX", ranges = IRanges(start = min(h9$BP), end= max(h9$BP)))

overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))

overlap$symbol = NA

for(i in 1:nrow(overlap)){
  overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
}


genes = overlap
genes = genes[-which(is.na(genes$symbol)),]
genes = genes[-which(genes$width < 100),]
#genes<-subset(genes,!duplicated(genes[,'name2']))
#genes<-genes[,c('txStart','txEnd','name2','strand')]
genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)


y0=y1=-1
for(i in 1:nrow(genes)){
  arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
         col='red')
  text(x=genes[i,'start'],y=-1,genes$symbol[i],cex=0.5, 
       col='red',srt=90)
  y0 = y0-0.015
  y1=y1-0.015
}

dev.off()




####loc10

  
  h10<- h[which(h$locus == 10),]
  
  
  cairo_pdf(file="~/Desktop/figs/sup1j.pdf", width = 10, height = 7)
  plot(h10$BP,-log10(h10$P.NI),type="p",
       col= 'red', cex=1.5,
       axes=T,ylab='-logP',xlab='Chromosome',ylim=c(-2,8), lwd=1)
  abline(h=7.30103)
  
  mapping = as.data.frame(org.Hs.egSYMBOL)
  
  gr = GRanges(seqnames = "chrX", ranges = IRanges(start = min(h10$BP), end= max(h10$BP)))
  
  overlap = as.data.frame(subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), gr))
  
  overlap$symbol = NA
  
  for(i in 1:nrow(overlap)){
    overlap$symbol[i] = mapping[which(mapping$gene_id == overlap$gene_id[i]),"symbol"]
  }
  
  
  genes = overlap
  genes = genes[-which(is.na(genes$symbol)),]
  genes = genes[-which(genes$width < 100),]
  #genes<-subset(genes,!duplicated(genes[,'name2']))
  #genes<-genes[,c('txStart','txEnd','name2','strand')]
  genes[,'strand2']<-ifelse(genes[,'strand']=='+',2,1)
  
  
  y0=y1=-1
  for(i in 1:nrow(genes)){
    arrows(x0=genes[i,'start'],y0=y0,x1=genes[i,'end'],y1=y1, code=genes[i,'strand2'], length=0.05, lwd=2, 
           col='red')
    text(x=genes[i,'start'],y=-1,genes$symbol[i],cex=0.5, 
         col='red',srt=90)
    y0 = y0-0.015
    y1=y1-0.015
  }
  
  dev.off()
  
  


