options(stringsAsFactors = FALSE)


#function to convert mouse to human gene names
#use mgi homolog list
#http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt (Mar-14-2020)

homology = read.table("./data/mgi_homologs.txt",sep = "\t", header = T)

convertMousetoHuman = function(x){
  human=c()
  for(i in 1:length(x)){
    id = homology[which((homology$Common.Organism.Name == "mouse, laboratory") & tolower(homology$Symbol) == tolower(x[i])),"HomoloGene.ID"]
    hum = homology[which((homology$Common.Organism.Name == "human") & homology$HomoloGene.ID == id),"Symbol"]
    human=append(human,hum)
  }
  human=unique(human)
  return(human)
 
}






#Get relevant genes from kda
kda = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv")
kda_m=read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv")
kda_f = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv")

nominal = kda[which(kda$hyper <= 0.05),"gene"]
nominal = append(nominal,kda_m[which(kda_m$hyper <= 0.05),"gene"])
nominal = append(nominal,kda_f[which(kda_f$hyper <= 0.05),"gene"])

nominal=unique(nominal) #1050 unique nominal key drivers

#convert mouse genes to human equivalent



nominal_hum = convertMousetoHuman(nominal) #904



##munge GTEx v7 morris output files

files = list.files("./results/flat/coloc/coloc_results_1mbp/")


all_coloc_results = data.frame(matrix(nrow=1,ncol=9))

for(f in files){
  #x = read.delim(paste0("~/Desktop/coloc_results_1mbp/",f),header = F,sep = "", stringsAsFactors = FALSE)
  x = read.delim(paste0("./results/flat/coloc/coloc_results_1mbp/",f),header = F,sep = "", stringsAsFactors = FALSE)
  x=x[-1,-1]
  dim(x) = c(9,length(x)/9)
  x = as.data.frame(x)
  x=t(x)
  colnames(x) = colnames(all_coloc_results)
  all_coloc_results = rbind(all_coloc_results,x)
}

all_coloc_results = all_coloc_results[-1,]

m = read.delim("./results/flat/coloc/morris_lead_snp_genes_1mbp.txt",header=T)
m2 = m[,c(4,5)]
m2 = unique(m2)

dupes = names(which(table(m2$ensembl_id) >1))

x = m2[which(m2$ensembl_id %in% dupes),]

all_coloc_results$gene = NA

all_coloc_results$gene = apply(all_coloc_results,1,function(x) m2[which(m2$ensembl_id == all_coloc_results$X3),"symbol"])
#m2 has duplicated ensembls

all_coloc_results = all_coloc_results[,-3]

colnames(all_coloc_results)= c("BMD","tissue","n","H0","H1","H2","H3","H4","gene")

coloc_over75 = all_coloc_results[which(as.numeric(all_coloc_results$H4) >= 0.75),]

coloc_gefos = read.delim("./results/flat/coloc/coloc_v7_all_results_int.txt")
genes = tolower(coloc_gefos$gene)
genes = append(genes, tolower(coloc_over75$gene))

###

coloc_gefos = coloc_gefos[,-1]
coloc_gefos = coloc_gefos[,-4]
coloc_over75 = coloc_over75[,c(1,2,9,3,4,5,6,7,8)]
colnames(coloc_gefos) = colnames(coloc_over75)
coloc = rbind(coloc_gefos, coloc_over75)
coloc = coloc[which(as.numeric(coloc$H4) >=0.75),]
coloc$gene = tolower(coloc$gene)

aa = coloc[order(coloc$gene, coloc$H4,decreasing = T), ] #sort by id and reverse of abs(value)
#aa = aa[ !duplicated(aa$gene), ]  
##
write.csv(aa,file = "./results/flat/coloc/all_coloc_greaterorover75",quote = F,row.names = F)


aa = read.csv("./results/flat/coloc/all_coloc_greaterorover75", stringsAsFactors = FALSE)
###
aa_nom = aa[which(aa$gene %in% tolower(nominal_hum)),]
#aa_ebmd = aa_nom[which(aa_nom$BMD == "BMD"),]
#aa_FN = aa_nom[which(aa_nom$BMD == "FNBMD"),]
#aa_LS = aa_nom[which(aa_nom$BMD == "LSBMD"),]

##
#annotate kda tables with coloc data
kda_all = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv")
kda_F = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv")
kda_M = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv")

#
kda_all$eBMD = kda_F$eBMD = kda_M$eBMD = NA
kda_all$FNBMD = kda_F$FNBMD = kda_M$FNBMD = NA
kda_all$LSBMD = kda_F$LSBMD = kda_M$LSBMD = NA
kda_all$PPH4 = kda_F$PPH4 = kda_M$PPH4 = NA

#kda_all
for(i in 1:nrow(kda_all)){
  if(tolower(kda_all$gene[i]) %in% aa_nom$gene & kda_all$hyper[i]<=0.05){
    x = subset(aa_nom, gene == tolower(kda_all$gene[i]))
    
    if("BMD" %in% x$BMD){kda_all$eBMD[i] = 1}
    if("FNBMD" %in% x$BMD){kda_all$FNBMD[i] = 1}
    if("LSBMD" %in% x$BMD){kda_all$LSBMD[i] = 1}
    kda_all$PPH4[i] = max(x$H4)
    
  }
}
########
#kda_F
#
for(i in 1:nrow(kda_F)){
  if(tolower(kda_F$gene[i]) %in% aa_nom$gene & kda_F$hyper[i]<=0.05){
    x = subset(aa_nom, gene == tolower(kda_F$gene[i]))
    
    if("BMD" %in% x$BMD){kda_F$eBMD[i] = 1}
    if("FNBMD" %in% x$BMD){kda_F$FNBMD[i] = 1}
    if("LSBMD" %in% x$BMD){kda_F$LSBMD[i] = 1}
    kda_F$PPH4[i] = max(x$H4)
    
  }
}


########
#kda_M
#
for(i in 1:nrow(kda_M)){
  if(tolower(kda_M$gene[i]) %in% aa_nom$gene & kda_M$hyper[i]<=0.05){
    x = subset(aa_nom, gene == tolower(kda_M$gene[i]))
    
    if("BMD" %in% x$BMD){kda_M$eBMD[i] = 1}
    if("FNBMD" %in% x$BMD){kda_M$FNBMD[i] = 1}
    if("LSBMD" %in% x$BMD){kda_M$LSBMD[i] = 1}
    kda_M$PPH4[i] = max(x$H4)
    
  }
}






#####ADD IMPC
kda_all$impc = kda_F$impc = kda_M$impc = NA


impc = read.csv("./results/flat/IMPC_BMD_Results.csv")
impc_gene = impc[which(impc$genotype_p_value<=0.05),"Gene"]

kda_all[which(tolower(kda_all$gene) %in% impc_gene),"impc"] = 1
kda_F[which(tolower(kda_F$gene) %in% impc_gene),"impc"] = 1
kda_M[which(tolower(kda_M$gene) %in% impc_gene),"impc"] = 1

write.csv(kda_all, file = "~/Desktop/kda_combined.csv", quote = F, row.names = F)
write.csv(kda_F, file = "~/Desktop/kda_F.csv", quote = F, row.names = F)

write.csv(kda_M, file = "~/Desktop/kda_M.csv", quote = F, row.names = F)


