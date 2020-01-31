#MAKE SURE TO MUS HUMAN GENES (COLOC)

##munge GTEx v7 morris output files

#files = list.files("~/Desktop/coloc_results_1mbp/")
files = list.files("./results/flat/coloc/coloc_results_1mbp/")

#tis_name = c()

#for (f in files){
#  tis_name = append(tis_name,strsplit(f,split = "coloc_results_")[[1]][2])
#}

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

#m = read.delim("~/Desktop/GWAS_project/morris_ukbb/morris_lead_snp_genes_1mbp.txt",header=T)
m = read.delim("./results/flat/coloc/morris_lead_snp_genes_1mbp.txt",header=T)


all_coloc_results$gene = NA
for(i in 1:nrow(all_coloc_results)){
  all_coloc_results$gene[i] = m[which(m$ensembl_id == all_coloc_results$X3[i]),"symbol"]
}
all_coloc_results = all_coloc_results[,-3]

colnames(all_coloc_results)= c("BMD","tissue","n","H0","H1","H2","H3","H4","gene")

coloc_over75 = all_coloc_results[which(as.numeric(all_coloc_results$H4) >= 0.75),]

#coloc_gefos = read.delim("~/Desktop/GWAS_project/project_extension/coloc_v7/coloc_v7_all_results_int.txt")
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
####
####
####
####
####
####
####
####
####
####
####

zhang = all_females_power8_BICOR


zhang$coloc_H0 = NA
zhang$coloc_H1 = NA
zhang$coloc_H2 = NA
zhang$coloc_H3 = NA
zhang$coloc_H4 = NA
zhang$coloc_tissue = NA
zhang$coloc_pheno = NA
zhang$coloc_FNBMD = 0
zhang$coloc_LSBMD = 0
zhang$coloc_eBMD = 0

for(i in which(tolower(aa$gene)  %in% tolower(zhang$gene))){
  gene = tolower(aa$gene[i])

  zhang[which(tolower(zhang$gene) == gene),"coloc_H0"] =  aa$H0[i]
  zhang[which(tolower(zhang$gene) == gene),"coloc_H1"] =  aa$H1[i]
  zhang[which(tolower(zhang$gene) == gene),"coloc_H2"] =  aa$H2[i]
  zhang[which(tolower(zhang$gene) == gene),"coloc_H3"] =  aa$H3[i]
  zhang[which(tolower(zhang$gene)== gene),"coloc_H4"] =  aa$H4[i]
  zhang[which(tolower(zhang$gene) == gene),"coloc_tissue"] = aa$tissue[i]
  zhang[which(tolower(zhang$gene) == gene),"coloc_pheno"] = aa$BMD[i]
}

coloc_FNBMD = coloc_gefos[which(coloc_gefos$BMD == "FNBMD"),]
coloc_LSBMD = coloc_gefos[which(coloc_gefos$BMD == "LSBMD"),]

zhang[which(tolower(zhang$gene) %in% tolower(coloc_FNBMD$gene)),"coloc_FNBMD"] = 1
zhang[which(tolower(zhang$gene) %in% tolower(coloc_LSBMD$gene)),"coloc_LSBMD"] = 1
zhang[which(tolower(zhang$gene) %in% tolower(coloc_over75$gene)),"coloc_eBMD"] = 1

colnames(zhang)[6] = "degree"

all_females_power8_BICOR = zhang


all_females_power8_BICOR$impc = 0
all_females_power8_BICOR[which(tolower(all_females_power8_BICOR$gene) %in% impc),"impc"] = 1





