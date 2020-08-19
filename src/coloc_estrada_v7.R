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


fnbmd = read.table("fnbmd_merged",stringsAsFactors = FALSE,as.is = TRUE, header=T)
lsbmd = read.table("lsbmd_merged",stringsAsFactors = FALSE,as.is = TRUE, header=T)

#get MAF from effect allele freq
idx = which(fnbmd$Freq.Allele1.HapMapCEU >0.5)
fnbmd$Freq.Allele1.HapMapCEU[idx] = 1-fnbmd$Freq.Allele1.HapMapCEU[idx]

idx = which(lsbmd$Freq.Allele1.HapMapCEU >0.5)
lsbmd$Freq.Allele1.HapMapCEU[idx] = 1-lsbmd$Freq.Allele1.HapMapCEU[idx]

gene_list = read.delim("estrada_lead_BAN_overlaps.txt",header = FALSE,as.is = TRUE,stringsAsFactors = FALSE, sep=" ")

colnames(gene_list) = c("rsid","hg19","chr","chrnum","pos","ensembl_id","symbol")

gene_list_sym = gene_list[,"symbol"]
gene_list_ens = gene_list[,"ensembl_id"]
gene_list_sym = unique(gene_list_sym)
gene_list_ens = unique(gene_list_ens)


##############

counter = 1
coloc_results = data.frame(matrix(NA, nrow = 40000, ncol = 9)) #empty table

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

  colnames(current_tissue_file) = c("ens_long", "snp_id_hg37", "tss_dist", "ma_samples", "ma_count","maf","pval", "slope", "slope_se", "lead_snp", "lead_snp_chrom", "lead_snp_pos","ensembl","gene_sym","lead_id19")



tissue_name=tissue_file

  analysis_genes = gene_list_ens[which(gene_list_ens %in% current_tissue_file$ensembl)] #genes that are in the tissue file
  #for gene in table
  
  for(gene in analysis_genes){
    gene_tissue_table = current_tissue_file[current_tissue_file$ensembl == gene,] #subset the tissue table into containing only the gene of interest

    fnbmd_cur = fnbmd[fnbmd$variant_id %in% gene_tissue_table$snp_id_hg37,]
    
    gtex_fn_cur = gene_tissue_table[gene_tissue_table$snp_id_hg37 %in% fnbmd_cur$variant_id,]
    
    if(length(which(duplicated(gtex_fn_cur$snp_id_hg37)))>0){
       gtex_fn_cur = gtex_fn_cur[-which(duplicated(gtex_fn_cur$snp_id_hg37)),]#some are duplicated because they are drawn under different lead snps
    }
    
    
    lsbmd_cur = lsbmd[lsbmd$variant_id %in% gene_tissue_table$snp_id_hg37,]
    
    gtex_ls_cur = gene_tissue_table[gene_tissue_table$snp_id_hg37 %in% lsbmd_cur$variant_id,]
    
    if(length(which(duplicated(gtex_ls_cur$snp_id_hg37)))>0){
       gtex_ls_cur = gtex_ls_cur[-which(duplicated(gtex_ls_cur$snp_id_hg37)),]#some are duplicated because they are drawn under different lead snps
     }

   	
    print(paste(counter,"of",length(analysis_genes),"genes",sep = " "))
    print(tissue_name)
    print(gene)
    tissue_n = as.numeric(tis[which(tis$Tissue == tissue_name),2]) #number of samples for the tissue
    
    gene.ls.coloc = list(pvalues=as.numeric(gtex_ls_cur$pval),N=as.numeric(tissue_n),type='quant',
                      snp=as.character(gtex_ls_cur$snp_id_hg37), MAF=as.numeric(gtex_ls_cur$maf))
    
    gene.fn.coloc = list(pvalues=as.numeric(gtex_fn_cur$pval),N=as.numeric(tissue_n),type='quant',
                       snp=as.character(gtex_fn_cur$snp_id_hg37), MAF=as.numeric(gtex_fn_cur$maf))

    fnbmd.coloc = list(pvalues=as.numeric(fnbmd_cur$P.value),N=32000,type='quant',
                       snp=as.character(fnbmd_cur$variant_id), MAF=as.numeric(fnbmd_cur$Freq.Allele1.HapMapCEU)) #
    
    lsbmd.coloc = list(pvalues=as.numeric(lsbmd_cur$P.value),N=32000,type='quant',
                        snp=as.character(lsbmd_cur$variant_id), MAF=as.numeric(lsbmd_cur$Freq.Allele1.HapMapCEU))
    
    
    coloc_FNBMD = try(coloc.abf(gene.fn.coloc,fnbmd.coloc))
    coloc_LSBMD = try(coloc.abf(gene.ls.coloc,lsbmd.coloc))



    if(class(coloc_FNBMD) != "try-error"){
      coloc_results[counter,] <<- c("FNBMD",tissue_name,gene,coloc_FNBMD$summary[1],coloc_FNBMD$summary[2],coloc_FNBMD$summary[3],
                                    coloc_FNBMD$summary[4],coloc_FNBMD$summary[5],coloc_FNBMD$summary[6])
    } else {coloc_results[counter,] <<- c("FNBMD",tissue_name,gene,NA,NA,NA,NA,NA,NA)}
    
    colnames(coloc_results) = c("pheno","tis_name","gene","nSNPs","H0","H1","H2","H3","H4")

    counter <<- counter + 1

  
  if(class(coloc_LSBMD) != "try-error"){
    coloc_results[counter,] <<- c("LSBMD",tissue_name,gene,coloc_LSBMD$summary[1],coloc_LSBMD$summary[2],coloc_LSBMD$summary[3],
                                  coloc_LSBMD$summary[4],coloc_LSBMD$summary[5],coloc_LSBMD$summary[6])
  } else {coloc_results[counter,] <<- c("LSBMD",tissue_name,gene,NA,NA,NA,NA,NA,NA)}
  
  colnames(coloc_results) = c("pheno","tis_name","gene","nSNPs","H0","H1","H2","H3","H4")
  
  counter <<- counter + 1
  
  }
}
  
  #return(coloc_results)

cur_tis = tis_list[slurm_arrayid,1]
print(cur_tis)
coloc_func(cur_tis)
colnames(coloc_results) = c("pheno","tis_name","gene","nSNPs","H0","H1","H2","H3","H4")
write.table(coloc_results,paste0("../coloc_results/","coloc_results_",cur_tis))

