library(IMPCdata)
library(PhenStat)

t1<-getIMPCDataset('WTSI','MGP_001','IMPC_DXA_001','IMPC_DXA_004_001','MGI:4364018')

t1<-subset(t1,t1$Genotype!="colony_id")

plot(as.vector(as.character(t1$Weight)),as.numeric(as.character(t1$Value)))

t1$Value<-as.numeric(as.character(t1$Value))
t1$Weight<-as.numeric(as.character(t1$Weight))

o1<-unique(t1$Genotype)
d2<-PhenList(t1, testGenotype=as.character(unique(t1$Genotype)[which(o1!='+/+')]),refGenotype='+/+',outputMessages = F)
d3<-testDataset(d2,depVariable='Value',method='MM',equation='withWeight', outputMessages = T,transformValues = T)

d3@analysisResults$model.output.summary
