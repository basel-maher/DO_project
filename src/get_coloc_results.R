#remove NAs and merge all coloc_all results

#morris eBMD GWAS data
#coloc_results_files =list.files(path= "~/Documents/projects//DO_project/results/flat/coloc/coloc_results_morris/", pattern = "coloc_results*")

coloc_results_files =list.files(path= "~/Documents/projects//DO_project/results/flat/coloc/coloc_morris_results_REV2/", pattern = "coloc_results*")
coloc_7_all_results = data.frame()
for (i in coloc_results_files){
  #result = read.delim(paste("~/Documents/projects//DO_project/results/flat/coloc/coloc_results_morris/",i,sep = ""),as.is = TRUE,stringsAsFactors = FALSE,header = FALSE,sep=" ")
  result = read.delim(paste("~/Documents/projects//DO_project/results/flat/coloc/coloc_morris_results_REV2/",i,sep = ""),as.is = TRUE,stringsAsFactors = FALSE,header = FALSE,sep=" ")
  result = result[-1,]
  result = result[which(is.na(result$V2)==FALSE),]
  coloc_7_all_results=rbind(coloc_7_all_results,result)
}

coloc_7_all_results = coloc_7_all_results[,-c(1)]
#gene_list = read.table("~/Documents/projects//DO_project/results/flat/coloc/morris_lead_BAN_overlaps.txt")
gene_list = read.table("~/Documents/projects//DO_project/results/flat/coloc/morris_lead_BAN_overlaps_REV2.txt")
gene_list = gene_list[,c(6,7)]
gene_list = unique(gene_list)


coloc_7_all_results$gene = apply(coloc_7_all_results,1,function(x) gene_list[which(gene_list[,1] == x[3]),2])                                     
                                     
colnames(coloc_7_all_results) = c("pheno","tissue","ensembl","nSNPs","H0","H1","H2","H3","H4","gene")
coloc_7_all_results$nSNPs = as.numeric(coloc_7_all_results$nSNPs)

coloc_7_all_results$gene = as.character(coloc_7_all_results$gene)
write.table(coloc_7_all_results,"~/Documents/projects/DO_project/results/flat/coloc/coloc_morris_v7_all_results_REV2.txt",sep = "\t",quote = FALSE)

coloc_7_all_results_over75 = coloc_7_all_results[which(coloc_7_all_results$H4 >=0.75),] 

morris_ovr75  = coloc_7_all_results_over75
write.table(coloc_7_all_results_over75,"~/Documents/projects/DO_project/results/flat/coloc/coloc_morris_v7_all_results_over75_REV2.txt",sep = "\t",quote = FALSE)


####repeat analysis for estrada GWAS


#coloc_results_files =list.files(path= "~/Documents/projects//DO_project/results/flat/coloc/coloc_results_estrada/", pattern = "coloc_results*")

coloc_results_files =list.files(path= "~/Documents/projects//DO_project/results/flat/coloc/coloc_estrada_results_REV2/", pattern = "coloc_results*")

coloc_7_all_results = data.frame()
for (i in coloc_results_files){
  result = read.delim(paste("~/Documents/projects//DO_project/results/flat/coloc/coloc_estrada_results_REV2/",i,sep = ""),as.is = TRUE,stringsAsFactors = FALSE,header = FALSE,sep=" ")
  result = result[-1,]
  result = result[which(is.na(result$V2)==FALSE),]
  coloc_7_all_results=rbind(coloc_7_all_results,result)
}

coloc_7_all_results = coloc_7_all_results[,-c(1)]
#gene_list = read.table("~/Documents/projects//DO_project/results/flat/coloc/estrada_lead_BAN_overlaps.txt")
gene_list = read.table("~/Documents/projects//DO_project/results/flat/coloc/estrada_lead_BAN_overlaps_REV2.txt")

#gene_list = gene_list[,c(6,7)]
gene_list = gene_list[,c(5,6)]
gene_list = unique(gene_list)


coloc_7_all_results$gene = apply(coloc_7_all_results,1,function(x) gene_list[which(gene_list[,1] == x[3]),2])                                     

colnames(coloc_7_all_results) = c("pheno","tissue","ensembl","nSNPs","H0","H1","H2","H3","H4","gene")
coloc_7_all_results$nSNPs = as.numeric(coloc_7_all_results$nSNPs)

coloc_7_all_results$gene = as.character(coloc_7_all_results$gene)
write.table(coloc_7_all_results,"~/Documents/projects/DO_project/results/flat/coloc/coloc_estrada_v7_all_results_REV2.txt",sep = "\t",quote = FALSE)

coloc_7_all_results_over75 = coloc_7_all_results[which(coloc_7_all_results$H4 >=0.75),] 

estrada_ovr75  = coloc_7_all_results_over75


fn = coloc_7_all_results_over75[which(coloc_7_all_results_over75$pheno=="FNBMD"),]
ls = coloc_7_all_results_over75[which(coloc_7_all_results_over75$pheno=="LSBMD"),]

write.table(fn,"~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_FNBMD_over75_REV2.txt",sep = "\t",quote = FALSE)
write.table(ls,"~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_LSBMD_over75_REV2.txt",sep = "\t",quote = FALSE)


coloc_ovr_75 = rbind(estrada_ovr75, morris_ovr75)
unique(coloc_ovr_75$gene)
