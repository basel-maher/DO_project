set.seed(8675309)
library(qtl2)
library(Mus.musculus)#seems to use mm10
library(topGO)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_REDO.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#
annot_file = read.csv("./results/flat/annot_file.csv", stringsAsFactors = FALSE)


###
#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file



#calculate ellipticity
a = cross_basic$pheno[,"ML"]
b = cross_basic$pheno[,"AP"]

e = sqrt((a^2 - b^2)/a^2)
e = log10(e)

e=b/a
#map ellipticity
ellipticity_scan = scan1(apr, e, k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])
e_peaks = find_peaks(ellipticity_scan, cross_basic$pmap, threshold=4, drop=1.5)






#eqtl correlations

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

for(i in 1:nrow(cors)){
  cors$name[i] = annot_file[which(annot_file$Gene.ID == cors$id[i]),"Gene.Name"]
}
cors$name = NA
cors$name = apply(cors, 1, function(x) annot_file[which(annot_file$Gene.ID == x[2]),"Gene.Name"])



top100 = unlist(cors[order(cors$abs, decreasing = T),"name"][2:101])#absolute top 100
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

