set.seed(8675309)
library(qtl2)
options(stringsAsFactors = F)



#####################
#Take all loci, get CI of locus +/- 250kb, get all eQTL within that locus, then compare top snps for phenos with top eqtl snps, define list of putatively causal genes

#qtl file
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

#remove FFP and soleus and MAT vol1_nonzero
#qtl_norm = qtl_norm[-c(1:3),]
#remove MAT_vol1_nonzero
#qtl_norm = qtl_norm[-(which(qtl_norm$lodcolumn == "MAT_VOL1_nonzero")),]

#define loci
qtl_norm$locus = 0

loc_idx = 1
qtl_loc = as.data.frame(matrix(nrow=nrow(qtl_norm), ncol=ncol(qtl_norm)))
colnames(qtl_loc) = colnames(qtl_norm)
for (i in c(1,2,3,4,8,10,16,"X")){
  sub = subset(qtl_norm, qtl_norm$chr == i)
  
  while(any(sub$locus == 0)){
    min_sub = min(sub$pos)
    idx = which((sub$pos >= min_sub) & (sub$pos <= min_sub + 1.5))
    sub[idx,"locus"] = loc_idx
    qtl_loc = rbind(qtl_loc, sub[idx,])
    sub = sub[-idx,]
    loc_idx = loc_idx + 1
  }
}

qtl_loc = qtl_loc[-which(is.na(qtl_loc)),]

#write.csv(qtl_loc, file = "./results/flat/qtl_loc", quote = FALSE,row.names = FALSE)

qtl_loc = read.csv("./results/flat/qtl_loc")


#for each locus, take genes that have eqtl , and are within locus CI +/- 250 kb
#any portion of gene start or end is within CI

#load eqtl
load("./results/Rdata/local_eqtl.Rdata")

eqtl_loc = as.data.frame(matrix(nrow=nrow(local_eqtl), ncol=ncol(local_eqtl)+1))
colnames(eqtl_loc) = colnames(local_eqtl)
colnames(eqtl_loc)[17] = "locus"

for(i in unique(qtl_loc$locus)){
  print(i)
  sub = subset(qtl_loc, qtl_loc$locus == i)
  min_loc = min(sub$ci_lo) - 0.25
  max_loc = max(sub$ci_hi) + 0.25
  
  sub_eqtl_1 = subset(local_eqtl, local_eqtl$chr == unique(sub$chr) & (local_eqtl$Start/1000000 >= min_loc) & local_eqtl$Start/1000000 <= max_loc)
  sub_eqtl_2 = subset(local_eqtl, local_eqtl$chr == unique(sub$chr) & (local_eqtl$End/1000000 >= min_loc) & local_eqtl$End/1000000 <= max_loc)
  
  sub_eqtl = merge(sub_eqtl_1, sub_eqtl_2)
  sub_eqtl = unique(sub_eqtl)
  sub_eqtl$locus = i
  
  eqtl_loc = rbind(eqtl_loc, sub_eqtl)
  
}
eqtl_loc = eqtl_loc[-which(is.na(eqtl_loc)),]



####
####
#load merge analysis objects

load("./results/Rdata/merge_top_local_eqtl.Rdata")
merge_top_eqtl = merge_top

load("./results/Rdata/merge_top_QTL.Rdata")
merge_top_qtl = merge_top

rm(merge_top)


#for each locus, take eqtl merge analyses in locus and colocalize with pheno merge analyses in locus


genes = c()
phenos_w_genes = c()
for(i in unique(qtl_loc$locus)){
  #print(i)
  sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
  sub_eqtl = subset(eqtl_loc, eqtl_loc$locus == i)
  
  phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
  
  eqtl_genes = paste0(sub_eqtl$lodcolumn, "_", unique(sub_eqtl$chr))
  
  
  for(j in phenos){
    for(k in eqtl_genes){
      print(k)
      if(any(merge_top_eqtl[[k]]$snp_id %in% merge_top_qtl[[j]]$snp_id)){
        print(k)
        genes = append(genes,k)
        phenos_w_genes = append(phenos_w_genes, j)
      }
    }
  }
}


unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
gene_names = unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
gene_names = unique(gene_names)

eqtl_loc[which(eqtl_loc$lodcolumn %in% genes),]


phenos_w_genes = as.data.frame(phenos_w_genes)
phenos_w_genes$gene = genes
phenos_w_genes$gene_name = NA

for(i in 1:nrow(phenos_w_genes)){
  phenos_w_genes$gene[i] = unlist(strsplit(phenos_w_genes$gene[i], "_"))[1]
  phenos_w_genes$gene_name[i] = eqtl_loc[which(eqtl_loc$lodcolumn == phenos_w_genes$gene[i]),"Gene.Name"]
}

#eqtl genes that are located within a phenotypic qtl and regulated by a local eqtl
phenos_w_genes




## number nonsynonymous variants in top qtl merge for a phenotype
nonsyn = list()
counter=1
for(i in unique(qtl_loc$locus)){
  #print(i)
  sub_qtl = subset(qtl_loc, qtl_loc$locus == i)

  phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
  
  
  for(j in phenos){
    counter=counter+1
    missense = (grep("missense",x=merge_top_qtl[[j]]$consequence))
    l = length(missense)
    z = merge_top_qtl[[j]]$snp_id[missense]
    print(paste0(j,":",l))
    nonsyn[[counter]] = z
    names(nonsyn)[counter] = j
    }
  }

nonsyn_frame = do.call(rbind, lapply(nonsyn, as.data.frame))




#SIFT: RSIDs uploaded to https://useast.ensembl.org/Tools/VEP


