options(stringsAsFactors = FALSE)
library(topGO)

#kda_analyses#
all = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4_REV.csv", stringsAsFactors = F)
all_f = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4_REV.csv", stringsAsFactors = F)
all_m = read.csv("./results/flat/key_driver_analysis_MALES_sft5_REV.csv", stringsAsFactors = F)

#bone gene superset
superset = read.delim("./results/flat/superset_in_networks.txt", stringsAsFactors = FALSE, header = FALSE)

superset = superset[,1]


#Analysis
genes = c(all[which(all$hyper<=0.05), "gene"], all_m[which(all_m$hyper<=0.05), "gene"], all_f[which(all_f$hyper<=0.05), "gene"])
genes = unique(genes)

length(which(tolower(genes) %in% tolower(superset)))
length(which(tolower(genes) %in% tolower(superset)==FALSE))

##are bone genes more highly connected in the network?
x = all #repeat for males and females
#x = all_m
#x = all_f

bone_df = x[which(tolower(x$gene) %in% tolower(superset)),]
not_bone_df = x[which(tolower(x$gene) %in% tolower(superset) == FALSE),]


#wilcox.test(bone_df$degree, not_bone_df$degree, alternative = "greater")#1.8e-4, 0.0283 males, 1.052e-5 females


wilcox.test(bone_df$num_neib, not_bone_df$num_neib, alternative = "greater")#2.621e-5, 0.005216 males, 5.542e-7 females


#we use this because the BAN was defined as being conencted to a bone gene with num_neib (3-step). 
###GO analysis
load("./results/Rdata/networks/geneModMemAnnot_power4.RData")
load("./results/Rdata/networks/geneModMemAnnot_m_power5.RData")
load("./results/Rdata/networks/geneModMemAnnot_f_power4.RData")

allGenes = unique(c(combat_annot$Gene.Name, combat_annot_f$Gene.Name, combat_annot_m$Gene.Name))

interesting.genes = not_bone_df$gene

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




