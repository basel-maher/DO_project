library(dplyr)
##merge analysis
##Do merge analysis for each of the significant QTL, in the widest confidence interval for each locus

#DOWNLOAD THE FILES. DOWNLOADED NOV. 14 2019
#download.file("https://ndownloader.figshare.com/files/18533342", "./data/CCdb/cc_variants.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609252", "./data/CCdb/mouse_genes_mgi.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609261", "./data/CCdb/mouse_genes.sqlite")
#

set.seed(8675309)
library(qtl2)

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

#load qtl mapping object
load("./results/Rdata/DO_qtl_scan_norm.Rdata")

annot_file = read.csv("./results/flat/annot_file.csv", stringsAsFactors = FALSE)



#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file

##
covar_eqtl = as.matrix(cross_eqtl$covar)
covar_eqtl[,"sex"] = (covar_eqtl[,"sex"] == "M")*1

covar_eqtl = covar_eqtl[,-1]#remove sac date as covar for now

covar_eqtl = apply(covar_eqtl,2,as.numeric) #make sure all cols are numeric
rownames(covar_eqtl) = rownames(cross_eqtl$covar)#make sure rownames match original cross file


##FIX QTL mapping for MAT nonzero vs full
#get qtl list, passed threshold
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

#define locus, 1Mbp
chroms = as.vector(unique(qtl_norm$chr)) #chroms
y <- c(1:19, "X","Y","MT")
chroms = chroms[order(match(chroms, y))]



qtl_out = data.frame(matrix(ncol=ncol(qtl_norm)+1))
locus = 1
#for a chromsome, subset qtl
for(i in chroms){
  df = subset(qtl_norm, qtl_norm$chr==i)
  df = df[order(df$pos),]
  
  
  df$locus = NA
  
  while(length(which(is.na(df$locus))) >=1){
    
    idx = which(is.na(df$locus))
    min_pos = min(df$pos[idx])
    df[idx,][which((df$pos[idx] - min_pos) <=1), "locus"] = locus
    locus = locus + 1
  }


  
  colnames(qtl_out) = colnames(df)
  qtl_out = rbind(qtl_out, df)
  
  for(j in 1:nrow(df)){
    x = df[j,]
    xx = df[which((abs(df$pos - df$pos[j])<=1) &(df$locus != df$locus[j]) & (df$pos!= df$pos[j])),]
    x = rbind(x,xx)
    
    if(nrow(x)>1){
      x$locus = locus
      qtl_out = rbind(qtl_out, x)
      locus = locus+1
    }
  }
}

qtl_out = qtl_out[-1,]

#remove singletons
rmv=c()

for(i in 1:nrow(qtl_out)){
  sub = subset(qtl_out, qtl_out$locus == qtl_out$locus[i])
  if(nrow(sub) == 1){
    x = which(qtl_out$chr == qtl_out$chr[i] & qtl_out$pos == qtl_out$pos[i] & qtl_out$lodcolumn == qtl_out$lodcolumn[i] )
    if(length(x) >1){
      print(length(x))
      rmv = c(rmv,i)
    }
  }
  

}
qtl_out = qtl_out[-rmv,]

#remove duplicated loci with different locus id

#which(duplicated(qtl_out[,1:8]))





#qtl_peaks = qtl_peaks_both_norm[which(qtl_peaks_both_norm$lod > qtl_peaks_both_norm$perm_thresh),]
#

merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")

for(i in unique(qtl_out$locus)){
  print(i)
  sub = subset(qtl_out, qtl_out$locus == i)
  chr = unique(sub$chr)
  start = min(sub$ci_lo)
  end = max(sub$ci_hi)
  
  variants_locus = query_variants(chr, start, end)
  genes_locus <- query_genes(chr, start, end)
  
  genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]
  
  if("pseudogene" %in% genes_locus$mgi_type){
    genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
  }
  
  if("miRNA gene" %in% genes_locus$mgi_type){
    genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
  }
  
  if("rRNA gene" %in% genes_locus$mgi_type){
    genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
  }
  
  
  merge[[i]] = list()
  
  for(j in unique(sub$lodcolumn)){
    print(j)
    out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,j], k_loco[[chr]],  addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],
                          query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
    
    out_name = paste0("./results/plots/merge_analysis/","chr_",chr,"_locus_",i,"_",j,".pdf")
    
    pdf(out_name) 
    
    plot_snpasso(out_snps$lod, out_snps$snpinfo, genes=genes_locus,drop_hilit=1.5)
    
    dev.off() 
    
    
    merge[[i]][[j]] = out_snps
    rm(out_snps)
    gc()
  }
}

plot_snpasso(merge[[1]]$uCT_Ct.TMD$lod, merge[[1]]$uCT_Ct.TMD$snpinfo, genes=genes_locus)

top <- top_snps(out_snps$lod, out_snps$snpinfo)

#out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno, k_loco[["1"]],  addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],
#                      query_func=query_variants,chr=1, start=154.26979, end=158.22213, keep_all_snps=TRUE)

plot_snpasso(out_snps$lod, out_snps$snpinfo)



plot_snpasso(out_snps$lod, out_snps$snpinfo, genes=genes_locus)

top <- top_snps(out_snps$lod, out_snps$snpinfo)
####
####
#plot bv/tv

x = as.data.frame(cross_basic$pheno)
x = x[order(x$uCT_BV.TV),]
x$row  = c(1:nrow(x))

p<-ggplot(data=x, aes(y=uCT_BV.TV, x=row, fill=row)) +geom_bar(stat="identity",width=0.5)+theme_minimal() +xlab("index (n=619)")+ylab("bone volume fraction (BV/TV, %)") + scale_fill_continuous() + theme(legend.position="none")

#bone volume fraction
merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")



###plot 
for(i in 1:nrow(qtl_out)){
  print(i)
  chr = qtl_out$chr[i]
  start = qtl_out$ci_lo[i]-2.5
  end = qtl_out$ci_hi[i]+2.5
  
  genes_locus <- query_genes(chr, start, end)
  
  genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]
  
  if("pseudogene" %in% genes_locus$mgi_type){
    genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
  }
  
  if("miRNA gene" %in% genes_locus$mgi_type){
    genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
  }
  
  if("rRNA gene" %in% genes_locus$mgi_type){
    genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
  }
  
  

  
 for(gene in unique(genes_locus$Name)){
    id = annot_file[which(annot_file$Gene.Name == gene),"Gene.ID"]
    if(isFALSE(length(id)==0)){
      
    
    
  
      if(id %in% colnames(cross_eqtl$pheno)){
        print(gene)
       
      
        minMarker = find_marker(cross_eqtl$pmap, chr = chr, pos =qtl_out$pos[i]-15 )
        maxMarker = find_marker(cross_eqtl$pmap, chr = chr, pos =qtl_out$pos[i]+15 )
      

      
        out_blup = scan1blup(apr[,chr],pheno = cross_eqtl$pheno[,id], kinship = k_loco[[chr]], addcovar =new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 2)
      
      
        out_name = paste0("./results/plots/gene_blups/","chr_",chr,"_pheno_",qtl_out$lodcolumn[i],"_",gene,".pdf")
      
        pdf(out_name) 
      
        i1 = which(rownames(out_blup) == minMarker)
        i2 = which(rownames(out_blup) == maxMarker)
      
        plot_coefCC(out_blup[c(i1:i2),], cross_eqtl$pmap,legend = "topleft",main = gene,scan1_output = subset(DO_qtl_scan_normal, lodcolumn=qtl_out$lodindex))
      
        dev.off() 
      
      
        rm(out_blup)
        gc()
    }
      
  }
}

}


###########DO MERGE ANALYSIS LIKE EQTL, BASED ON CI ONLY.##########

merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")

for(i in 1:nrow(qtl_norm)){
  print(i)
  merge[[i]] = list()
  pheno = qtl_norm$lodcolumn[i]
  chr = qtl_norm$chr[i]
  start = qtl_norm$ci_lo[i]
  end = qtl_norm$ci_hi[i]
  #use same covars as qtl mapping.sex,age,BW and gen
  out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,pheno], k_loco[[chr]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                        query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
  
}




#for each merge analysis, take the snps that are within 15% LODs of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)



#merge_top_df = bind_rows(merge_top, .id = "column_label")



load("~/Desktop/merge_analysis/merge_top_local_eqtl.Rdata")
merge_top_eqtl = merge_top

summary(unlist(lapply(merge_top, function(x) nrow(x))))

length(which(unlist(lapply(merge_top, function(x) nrow(x))) == 1))

load("~/Desktop/merge_analysis/merge_top_QTL.Rdata")

summary(unlist(lapply(merge_top, function(x) nrow(x))))

which(unlist(lapply(merge_top, function(x) nrow(x))) == 1)
merge_top[[24]]

##
mmarg = maxmarg(pr, map=cross_basic$pmap, chr=16, pos=22.89766, minprob=0.51, return_char = T)
#m_full = maxmarg(pr, minprob=0.5)

covar = cbind(covar, as.data.frame(mmarg))

mmarg_model_matrix = model.matrix(~factor(covar$mmarg))[,-1]

covar_mmarg = merge(covar, mmarg_model_matrix, by="row.names", all = TRUE)
rownames(covar_mmarg) = covar_mmarg$Row.names

plot_pxg(mmarg, log10(cross_basic$pheno[,"bending_frax_load"]),sort = T, SEmult = 1)

chr=16
start = 22.31704
end = 23.42906

out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,"uCT_pMOI"], k_loco[["X"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                      query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)

#UNCHS048601
plot_scan1(DO_qtl_scan_normal, map = cross_basic$pmap["X"][800:1700], lodcolumn = "uCT_pMOI",chr = "X")

out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,"bending_max_load"], k_loco[[16]],  addcovar =  covar_snp[,c(2,3,4,7:16,18:20)],Xcovar=Xcovar,
                      query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


out_blup = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_mmarg[,c(2,3,4,7:16,18:52)],cores = 2)

plot_coefCC(out_blup, cross_basic$pmap,scan1_output = subset(DO_qtl_scan_normal, lodcolumn="bending_max_load"))


top = top_snps(out_snps$lod, out_snps$snpinfo, drop=1.5)
top[order(top$lod,decreasing = T),]

#cast private 23.22667
mmarg = maxmarg(pr, map=cross_basic$pmap, chr=16, pos=23.22667, minprob=0.51, return_char = T)
mmarg_model_matrix = model.matrix(~factor(covar$mmarg))[,-1]

covar_mmarg = merge(covar, mmarg_model_matrix, by="row.names", all = TRUE)
rownames(covar_mmarg) = covar_mmarg$Row.names


out_blup_bland = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar[,c(1,2,3,6:15)],cores = 2)
max_scan1(out_blup_bland, map = cross_basic$pmap, lodcolumn = NULL)


#cast private
snpinfo <- data.frame(chr=c("16"),
                      pos=c(23.22667),
                      sdp=32,
                      snp=c("rs241887612"), stringsAsFactors=FALSE)

#cast private right over qtl peak
snpinfo <- data.frame(chr=c("16"),
                      pos=c(23.26738),
                      sdp=32,
                      snp=c("rs217130452"), stringsAsFactors=FALSE)
#rs4165304

snpinfo <- data.frame(chr=c("16"),
                      pos=c(22.89766),
                      sdp=36,
                      snp=c("rs4165304"), stringsAsFactors=FALSE)

#no snp near peak 


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`16`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
if(covar_snp$A.rs217130452[i] >0.6){
  covar_snp$alleleA[i] = 1
} else{if(covar_snp$B.rs217130452[i] >0.6){
  covar_snp$alleleA[i] = 0
}}
}

out_blup_rs4165304_factor = scan1coef(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_snp[,c(2,3,4,7:16,19)],cores = 2)
out_blup_rs217130452_factor = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_snp[,c(2,3,4,7:16,19)],cores = 2)
out_blup_rs241887612_factor = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_snp[,c(2,3,4,7:16,19)],cores = 2)


##
out_blup_rs4165304 = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_snp[,c(2,3,4,7:19)],cores = 2)

out_blup_rs217130452 = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_snp[,c(2,3,4,7:19)],cores = 2)

out_blup_rs241887612 = scan1blup(apr[,16],pheno = cross_basic$pheno[,"bending_max_load"], kinship = k_loco[[16]], addcovar =  covar_snp[,c(2,3,4,7:19)],cores = 2)

minMarker = find_marker(cross_basic$pmap, chr = 16, pos =23.26738-10 )
maxMarker = find_marker(cross_basic$pmap, chr = 16, pos =23.26738+10 )

i1 = which(rownames(out_blup) == minMarker)
i2 = which(rownames(out_blup) == maxMarker)

plot_coefCC(out_blup_rs241887612_factor[c(i1:i2),], cross_basic$pmap,scan1_output = subset(DO_qtl_scan_normal, lodcolumn="bending_max_load"))


out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,"bending_max_load"], k_loco[[16]],  addcovar =  covar_snp[,c(2,3,4,7:19)],Xcovar=Xcovar,
                      query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


##

genes_locus <- query_genes(chr, start, end)

plot_snpasso(out_snps$lod, out_snps$snpinfo,genes = genes_locus)


##eqtl
out_snps_eqtl <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,"MSTRG.9880"], k_loco[[chr]],  addcovar = covar_eqtl[,c(1,10:57)],Xcovar=Xcovar,
                      query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)

top_snps(out_snps_eqtl$lod, out_snps_eqtl$snpinfo, drop=1.5)

plot_snpasso(out_snps_eqtl$lod, out_snps_eqtl$snpinfo,genes = genes_locus)

plot_pxg(mmarg, cross_eqtl$pheno[,"MSTRG.9880"],sort = T)

x= scan1(pr, cross_eqtl$pheno[,"MSTRG.9880"],kinship = k_loco[[chr]])
plot_scan1(x,cross_eqtl$pmap,chr="2")




for(i in 1:length(merge_top)){
  if("rs4165304" %in% merge_top[[i]]$snp_id){
    print(names(merge_top)[i])
  }
}

"rs4165304" %in% merge_top[["ENSMUSG00000000247_2"]]$snp_id

m = find_marker(cross_eqtl$pmap, chr = chr, pos = 22.89766)




"rs4165304" %in% merge_top[["ENSMUSG00000000247_2"]]$snp_id
x=load("~/Desktop/merge_top_local_eqtl_0.3drop.Rdata")


#####################
#Take all loci, get CI of locus +/- 250kb, get all eQTL within that locus, then compare top snps for phenos with top eqtl snps, define list of putatively causal genes

#qtl file
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

#remove FFP and soleus
qtl_norm = qtl_norm[-c(1:3),]
#remove MAT_vol1_nonzero
qtl_norm = qtl_norm[-(which(qtl_norm$lodcolumn == "MAT_VOL1_nonzero")),]

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

write.csv(qtl_loc, file = "./results/flat/qtl_loc", quote = FALSE,row.names = FALSE)




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

load("~/Desktop/merge_analysis/merge_top_local_eqtl.Rdata")
merge_top_eqtl = merge_top

load("~/Desktop/merge_analysis/merge_top_QTL.Rdata")
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


#herit
hsq = est_herit(cross_basic$pheno, kinship = k, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")])

