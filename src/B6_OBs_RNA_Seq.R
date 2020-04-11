#created by Charles Farber


# read in all files

files<-list.files('./data/GSE54461_OB_RNA_SEQ_RAW')[-c(28)] #use brackets
# to exclude non-RNA-seq data files in the same folder

all.ob.data<-read.delim(paste0('./data/GSE54461_OB_RNA_SEQ_RAW/',files[1]),header=T,sep='\t')
head(all.ob.data)
all.ob.data1<-all.ob.data[,c(1,7)]
head(all.ob.data1)

# read in files and created aggregated data file
for(i in 2:length(files)){
  all.ob.data1[,(i+1)]<-read.delim(paste0('./data/GSE54461_OB_RNA_SEQ_RAW/',files[i]),header=T, sep='\t')[,7]
  print(i)
}

# label columns
all.ob.data1[1:5,]
days<-c('d2','d2','d2','d4','d4','d4','d6','d6','d6','d8','d8','d8','d10','d10','d10','d12','d12','d12','d14','d14','d14','d16','d16','d16','d18','d18','d18')
colnames(all.ob.data1)<-c('Gene',days)

# read in ensembl gene annotation
library(biomaRt)
mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",host='www.ensembl.org')
listAttributes(mouse)[1:100,]

gene.ann<-getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", values = as.character(all.ob.data[,1]),mart = mouse)

gene.ann<-read.csv('./data/GSE54461_OB_RNA_SEQ_RAW/mart_export.txt',header=T)
gene.ann<-gene.ann[,-2]
head(gene.ann);dim(gene.ann)
gene.ann<-subset(gene.ann, !duplicated(gene.ann[,1]))
#genes<-c('Bicc1')
  # select gene or genes
#grep('ENSMUST00000182158',gene.ann[,1])
#grep('A430106G13Rik',gene.ann[,2])
#dat1<-c()
genes<-c('Qsox1')#A930004D18Rik
genes<-unique(genes)
rm(dat1)
dat1<-c()
for(i in 1:length(genes)){
  t1<-gene.ann[which(gene.ann[,2]==genes[i]),1]
  t2<-which(as.character(all.ob.data[,1])==as.character(t1))
  dat1<-c(dat1,t2)
}
dat1
# average replicates
temp1<-all.ob.data1[dat1,]
temp1$Gene
class(temp1$variable)
#rownames(temp1)<-temp1[,1]
#temp1<-temp1[,-1]
#t(temp1)
temp1<-t(temp1)
colnames(temp1)<-temp1[1,]
temp1<-temp1[-1,]
temp1<-data.frame(temp1)
temp1
temp1$day<-c(rep('d2',3),rep('d4',3),rep('d6',3),rep('d8',3),rep('d10',3),rep('d12',3),
             rep('d14',3),rep('d16',3),rep('d18',3))

temp1$day<-factor(temp1$day,levels=c('d2','d4','d6','d8','d10','d12','d14','d16','d18'))
library(reshape2)
temp2<-melt(temp1,id.vars='day');head(temp2)
class(temp2$variable)
library(ggplot2)
library(plyr)

temp3<-temp2[,c(2,1,3)]
colnames(temp3)<-c('Gene','Day','value')
temp3$value<-log2(as.numeric(as.character(temp3$value))+1)

temp4<-ddply(temp3, c('Gene','Day'), summarise,
      mean = mean(value), sd = sd(value),
      sem = sd(value)/sqrt(length(value)))


temp4$newGene<-mapvalues(temp4$Gene, from=unique(temp4$Gene), to=toupper(genes))

save(temp4,file = "./results/Rdata/rasd1_calvarial.Rdata")

cairo_pdf(filename="~/Desktop/figs/6A.pdf", width = 10, height = 7)
ggplot(temp4, aes(x=Day, y=mean,ymin=mean-temp4$sem, ymax=mean+temp4$sem, colour=newGene,Group=newGene)) + 
  geom_errorbar(aes(ymin=mean-temp4$sem, ymax=mean+temp4$sem), width=0.5, size=1.25) +
  geom_line(aes(group=newGene), size=0.5) +
  geom_point() + scale_color_discrete(name="Gene",
                                      breaks=as.character(unique(temp4$Gene)),
                                      labels=genes,l=60) +
  scale_y_continuous() + facet_wrap(~newGene, scales='free_y') 
dev.off()
