library(coloc)
options(stringsAsFactors=F)

# grab the array id value from the environment variable passed from sbatch

args <- commandArgs(trailingOnly=TRUE)
print(args)
slurm_arrayid = as.numeric(args[1])
print(slurm_arrayid)

##########


tis_list = read.delim("tissue_list_v7", header = F, stringsAsFactors = F)
tissue_sample_number = read.csv("GTEx_summary_v7.csv",stringsAsFactors = FALSE,as.is = TRUE)
renamer = function(x){
  x[1] = gsub("- ","",x[1])#replace dash followed by space with nothign
  x[1] = gsub(" ","_",x[1])#replace space with underscore
  x[1] = gsub("\\(","",x[1])#remove parentheses
  x[1] = gsub("\\)","",x[1])
  return(x)
}

tis = apply(tissue_sample_number,1,renamer)
tis = t(tis)
tis = as.data.frame(tis,stringsAsFactors = FALSE)#tis is tissue_sample_number but with names changed


bmd = read.delim("./Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE,as.is = TRUE, sep = "\t")
bmd$id19_new = paste0(gsub(x = bmd$SNPID, pattern = ":", replacement = "_"),"_","b37")


#get MAF from effect allele freq
idx = which(bmd$EAF >0.5)
bmd$EAF[idx] = 1-bmd$EAF[idx]


gene_list = read.delim("morris_lead_BAN_overlaps.txt",header = FALSE,as.is = TRUE,stringsAsFactors = FALSE, sep=" ")

colnames(gene_list) = c("hg19","hg19_new","rsid","chr","pos","ensembl_id","symbol")

gene_list_sym = gene_list[,"symbol"]
gene_list_ens = gene_list[,"ensembl_id"]
gene_list_sym = unique(gene_list_sym)
gene_list_ens = unique(gene_list_ens)


##############

counter = 1
coloc_results = data.frame(matrix(NA, nrow = 20000, ncol = 9)) #empty table

coloc_func = function(tissue_file){
  current_tissue_file = read.delim(paste0("../merged_eqtl_200k/",tissue_file),stringsAsFactors = FALSE,header = FALSE, sep=" ")#readtissue file


temp <- gsub(pattern="\t", replacement=" ", current_tissue_file[,1], fixed = TRUE)

# strsplit, make dataframe
temp <- do.call(rbind.data.frame, strsplit(temp, split = " "))

# output, merge columns
output <- cbind(temp,
                current_tissue_file[, c(2:ncol(current_tissue_file))])
 
current_tissue_file = output
rm(output)

  colnames(current_tissue_file) = c("ens_long", "snp_id_hg37", "tss_dist", "ma_samples", "ma_count","maf","pval", "slope", "slope_se", "lead_snp", "lead_snp_chrom", "lead_snp_pos","ensembl","gene_sym","lead_id19","lead_id19_new")



tissue_name=tissue_file

  analysis_genes = gene_list_ens[which(gene_list_ens %in% current_tissue_file$ensembl)] #genes that are in the tissue file
  #for gene in table
  
  for(gene in analysis_genes){
    gene_tissue_table = current_tissue_file[current_tissue_file$ensembl == gene,] #subset the tissue table into containing only the gene of interest

    bmd_cur = bmd[bmd$id19_new %in% gene_tissue_table$snp_id_hg37,]
    gtex_cur = gene_tissue_table[gene_tissue_table$snp_id_hg37 %in% bmd_cur$id19_new,]
    
    if(length(which(duplicated(gtex_cur$snp_id_hg37)))>0){
      gtex_cur = gtex_cur[-which(duplicated(gtex_cur$snp_id_hg37)),]#some are duplicated because they are drawn under different lead snps
    }

   	
    print(paste(counter,"of",length(analysis_genes),"genes",sep = " "))
    print(tissue_name)
    print(gene)
    tissue_n = as.numeric(tis[which(tis$Tissue == tissue_name),2]) #number of samples for the tissue
    
    gene.coloc = list(pvalues=as.numeric(gtex_cur$pval),N=as.numeric(tissue_n),type='quant',
                      snp=as.character(gtex_cur$snp_id_hg37), MAF=as.numeric(gtex_cur$maf))
    
    bmd.coloc = list(pvalues=as.numeric(bmd_cur$P.NI),N=bmd_cur$N[1],type='quant',
                       snp=as.character(bmd_cur$id19_new), MAF=as.numeric(bmd_cur$EAF)) #
    
    
    coloc_BMD = try(coloc.abf(gene.coloc,bmd.coloc))
    
    print(nrow(gene.coloc))
    print(nrow(bmd.coloc))

    print(class(coloc_BMD))

    if(class(coloc_BMD) != "try-error"){
      coloc_results[counter,] <<- c("BMD",tissue_name,gene,coloc_BMD$summary[1],coloc_BMD$summary[2],coloc_BMD$summary[3],
                                    coloc_BMD$summary[4],coloc_BMD$summary[5],coloc_BMD$summary[6])
    } else {coloc_results[counter,] <<- c("BMD",tissue_name,gene,NA,NA,NA,NA,NA,NA)}
    
    colnames(coloc_results) = c("pheno","tis_name","gene","nSNPs","H0","H1","H2","H3","H4")

    counter <<- counter + 1

  }
  
  #return(coloc_results)
}
cur_tis = tis_list[slurm_arrayid,1]
print(cur_tis)
coloc_func(cur_tis)
colnames(coloc_results) = c("pheno","tis_name","gene","nSNPs","H0","H1","H2","H3","H4")
write.table(coloc_results,paste0("../coloc_BAN_results/","coloc_results_",cur_tis))

