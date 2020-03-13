#eqtl correlations wwith qsox1
library(topGO)

#annot file
annot_file = read.csv("~/Documents/projects/DO_project/results/flat/annot_file.csv")

#load eqtl cross file
load(file = "./results/Rdata/cross_eqtl_REDO.Rdata")





cors = c()
#get all correlations with Qsox1, top 100, GO analysis
for(i in 1:ncol(cross_eqtl$pheno)){
  print(i)
  x = cor(cross_eqtl$pheno[,"MSTRG.1311"], cross_eqtl$pheno[,i],use = "p", method = "p")
  cors = append(cors,x)
}
names(cors) = colnames(cross_eqtl$pheno)
x = sort(abs(cors),decreasing = T)[1:101]
n = names(x)

cors = as.data.frame(cors)
cors$id = rownames(cors)
cors$abs = abs(cors$cors)

cors$name = NA
cors$name = apply(cors, 1, function(x) unlist(annot_file[which(annot_file$Gene.ID == x[2]),"Gene.Name"]))

top100 = unlist(cors[order(cors$cors, decreasing = T),"name"][2:101])#pos cor top 100


##topGO
allGenes = unname(unlist(cors$name))
interesting.genes = unname(top100)

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

