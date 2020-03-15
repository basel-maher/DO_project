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

#files = list.files("~/Desktop/coloc_results_1mbp/")
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


all_coloc_results$gene = NA
for(i in 1:nrow(all_coloc_results)){
  all_coloc_results$gene[i] = m[which(m$ensembl_id == all_coloc_results$X3[i]),"symbol"]
}
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
aa = aa[ !duplicated(aa$gene), ]  
##
#write.csv(aa,file = "./results/flat/coloc/all_coloc_greaterorover75",quote = F,row.names = F)


aa = read.csv("./results/flat/coloc/all_coloc_greaterorover75", stringsAsFactors = FALSE)
###
nominal_hum[which(tolower(nominal_hum) %in% aa$gene)]







