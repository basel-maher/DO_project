#figures/tables for paper

set.seed(8675309)
library(qtl2)
library(ggplot2)
library(gridExtra)
library(grid)
library(patchwork)
library(gtable)
library(ggrepel)
library(ggplotify)
#
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

##########
##########

#1B)
#freq of founder genotype per individual
af_ind <- calc_geno_freq(apr, "individual")#marker

mean_allele_freq = as.data.frame(matrix(nrow=8, ncol=1))

for(i in 1:8){
  mean_allele_freq[i,] = mean(af_ind[,i])
}
mean_allele_freq$allele = c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")
mean_allele_freq$x="A"

p=ggplot(mean_allele_freq, aes(fill=allele, y=V1, x=allele)) + 
  geom_bar(position="stack", stat="identity")#position=fill for percentage in stack

#x=x to stack, x=allele for normal bar plot

p+ scale_fill_manual(values=as.vector(CCcolors))

# 
# par(mar=c(4.1, 4.1, 0.6, 0.6))
# 
# par(mfrow=c(4,2))
# for(i in 1:8){ hist(af_ind[,i], main=NULL, breaks=30,xlab = "geno freq by ind.",cex.lab=1.25, cex.axis=1.25)
#   abline(v=mean(af_ind[,i]),col="red",lwd=3)
# }



#1C)
#boxplot BV/TV
#make dataframe with bv/tv and sex
df_bvtv = as.data.frame(cross_basic$pheno[,"uCT_BV.TV"])
df_bvtv$pheno = "BV/TV"
df_bvtv$ind = rownames(df_bvtv)
df_bvtv$M = covar[,"sex"]
colnames(df_bvtv)[1] = "val"
df_2 = df_bvtv
df_bvtv$pheno = paste0(df_bvtv$pheno,"_",df_bvtv$M)
#M=1

df = rbind(df_bvtv,df_2)


p1 <- ggplot(df, aes(x=pheno, y=val)) + 
  geom_boxplot()



df_str = as.data.frame(cross_basic$pheno[,"bending_max_load"])
df_str$pheno = "max_load"
df_str$ind = rownames(df_str)
df_str$M = covar[,"sex"]
colnames(df_str)[1] = "val"

df_2 = df_str
df_str$pheno = paste0(df_str$pheno,"_",df_str$M)
#M=1

df = rbind(df_bvtv,df_str)


p2 <- ggplot(df, aes(x=pheno, y=val)) + 
  geom_boxplot()

p1 | p2

df = rbind(df_bvtv, df_str)
#df_2 = df

#df$pheno = paste0(df$pheno,"_",df$M)
#df = rbind(df,df_2)

p2 <- ggplot(df, aes(x=pheno, y=val)) + 
  geom_boxplot()
p2
# 
# p
#1D) heritability

norm_pheno = as.data.frame(log10(cross_basic$pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,52,54:58,61:70,72,74,76)]))
pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,36,42,50,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}


new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.


#heritability
h = est_herit(pheno = pheno_combined, kinship = k,addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")] )

h_df = as.data.frame(h)
h_df$pheno = rownames(h_df)
# #males
# m = which(covar[,"sex"] == 1)
# hsq_m = est_herit(pheno = cross_basic$pheno[m,], kinship = k, addcovar = covar[,c("age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")], max_batch=1)
# hsq_m_df = as.data.frame(hsq_m)
# 
# #females
# f = which(covar[,"sex"] == 0)
# hsq_f = est_herit(pheno = cross_basic$pheno[f,], kinship = k, addcovar = covar[,c("age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")], max_batch=1)
# hsq_f_df = as.data.frame(hsq_f)


#hsq_table= merge(hsq_df, hsq_m_df, by=0)
#rownames(hsq_table) = hsq_table$Row.names
#hsq_table = hsq_table[,-1]

#hsq_table= merge(hsq_table, hsq_f_df, by=0)

#write.table(hsq_table, file="~/Desktop/hsq.csv", quote=F, sep=",")


#grid.table(hsq_table)
hsq_table = h_df[-which(h_df$pheno %in% c("glucose", "soleus_weight","adiposity","RFP","GFP","BFP","FFP","gastroc_weight","MAT_VOL1_nonzero","MAT_VOL2_nonzero","MAT_VOL3_nonzero","MAT_VOL4_nonzero")),]

hsq_table = hsq_table[order(hsq_table$h),]

hsq_table$pheno <- factor(hsq_table$pheno, levels = hsq_table$pheno)

p3<-ggplot(hsq_table, aes(x=pheno, y=h)) + 
  geom_point() + coord_flip()
p3



#wrap_plots(p1,p,p2,p3)
#patchwork+plot_annotation(tag_levels = "A")

# x = as.data.frame(hsq)
# x$pheno = rownames(x)
# x$sex = "B"
# 
# y = hsq_m_df
# y$pheno = rownames(y)
# y$sex = "M"
# colnames(y)[1] = "hsq"
# 
# z = hsq_f_df
# z$pheno = rownames(z)
# z$sex = "F"
# colnames(z)[1] = "hsq"
# 
# x = rbind(x,y)
# x = rbind(x,z)
# x = x[-which(x$pheno %in% c("sex", "sac_time", "age_at_sac_days","Mouse.ID")),]
# 
# 
# p<-ggplot(x, aes(x=pheno, y=hsq)) + 
#   geom_point(aes(color=sex, shape=sex)) + coord_flip()
#####
#####
#####
p1=grid::textGrob("PLACEHOLDER FOR A")
p1 = wrap_elements(p1)

patchwork = (p1|p|p2)/p3

layout = "
AABC
DDDD
"

p1+p+p2+p3+plot_layout(design = layout)+plot_annotation(tag_levels = "A")



#
#
#
#

##Figure 2
#2A
#scatterplot with traits, across all chroms, color traits differently, LOD score, labels
#load in significant qtl
qtl = read.csv("./results/flat/qtl_loc",stringsAsFactors = F)
qtl$pheno_name = c("ML","Ma.Ar","Tt.Ar","TMD","Ct.Por","pMOI","Imax","Ct.Ar/Tt.Ar","MAT_VOL1","ML","Ma.Ar","Ma.Ar","Tt.Ar","Ct.Ar/Tt.Ar","BMD","Disp. @ frax","Disp. @ max load","Tot.Work","WPY","TMD","Max Load","Frax. Load","Ct.Ar","pMOI","Imax","Imin","Tb.Sp","Tb.N","Ct.Th")
qtl[which(qtl$chr=="X"),"chr"]="20"
qtl$chr = as.numeric(qtl$chr)
chr_lengths = as.data.frame(matrix(nrow=2,ncol=1))
chr_lengths[1,1]=0
chr_lengths[2,1]=195
chr_lengths[3,1]=chr_lengths[2,1]+182
chr_lengths[4,1]=chr_lengths[3,1]+160
chr_lengths[5,1]=chr_lengths[4,1]+157
chr_lengths[6,1]=chr_lengths[5,1]+152
chr_lengths[7,1]=chr_lengths[6,1]+150
chr_lengths[8,1]=chr_lengths[7,1]+145
chr_lengths[9,1]=chr_lengths[8,1]+129
chr_lengths[10,1]=chr_lengths[9,1]+125
chr_lengths[11,1]=chr_lengths[10,1]+131
chr_lengths[12,1]=chr_lengths[11,1]+122
chr_lengths[13,1]=chr_lengths[12,1]+120
chr_lengths[14,1]=chr_lengths[13,1]+120
chr_lengths[15,1]=chr_lengths[14,1]+125
chr_lengths[16,1]=chr_lengths[15,1]+104
chr_lengths[17,1]=chr_lengths[16,1]+91
chr_lengths[18,1]=chr_lengths[17,1]+95
chr_lengths[19,1]=chr_lengths[18,1]+91
chr_lengths[20,1]=chr_lengths[19,1]+61
chr_lengths[21,1]=chr_lengths[20,1]+171

for(i in 1:nrow(qtl)){
  qtl$pos2[i] = chr_lengths[qtl$chr[i],] + qtl$pos[i]
}

qtl$pos2 = as.numeric(qtl$pos2)

#Breaks for background rectangles
rects <- data.frame(start = c(chr_lengths[1:20,1]), end = c(chr_lengths[2:21,1]), col =c(rep(c(1:2),10)))


ggplot() + 
  geom_rect(data = rects, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = col),fill=rep(c("grey","white"),10), alpha = 0.4) +
  geom_point(data = qtl, aes(x=pos2,y=lod))+ labs(x = "Chromosome", y = "LOD") + 
  scale_x_continuous(breaks=c(195/2,
                               (195+377)/2,
                               (377+537)/2,
                               (537+694)/2,
                               (694+846)/2,
                               (846+996)/2,
                               (996+1141)/2,
                               (1141+1270)/2,
                               (1270+1395)/2,
                               (1395+1528)/2,
                               (1526+1648)/2,
                               (1648+1768)/2,
                               (1768+1888)/2,
                               (1888+2013)/2,
                               (2013+2117)/2,
                               (2117+2208)/2,
                               (2208+2303)/2,
                               (2303+2394)/2,
                               (2394+2455)/2,
                               (2455+2626)/2),
                     limits=c(0,2626),labels = c(1:19,"X"),expand = c(0,0))+
  geom_text_repel(data=qtl,aes(x=pos2,y=lod,label=pheno_name),size = 3.5)+
  theme(legend.position = "none",panel.grid = element_blank())
#color by group?
  #geom_label_repel


###2B
#chr3 Ma.Ar QTL Trace

load("./results/Rdata/DO_qtl_scan_norm.Rdata")

maar3_blup = scan1blup(apr[,3], pheno_combined[,"uCT_Ma.Ar"], k_loco[["3"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
#save(maar3_blup,file = "./results/Rdata/maar_chr3_blup.Rdata")

plot_coefCC(maar3_blup[1000:4000,], cross_basic$pmap,legend = "topright",scan1_output = subset(DO_qtl_scan_normal, lodcolumn="uCT_Ma.Ar"), main = "Ma.Ar", legend_ncol=1,top_panel_prop = 0.6)
# 




##2C

query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")



start = 66.56421 - 0.250
end = 69.96983 + 0.250
chr = 3
out_snps_maar_3 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ma.Ar"], k_loco[["3"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                           query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=T)


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

top_maar3 <- top_snps(out_snps_maar_3$lod, out_snps_maar_3$snpinfo, drop = 0.15 * max(out_snps_maar_3$lod))
top_maar3[order(top_maar3$lod,decreasing = T),]


plot_snpasso(out_snps_maar_3$lod, out_snps_maar_3$snpinfo, genes = genes_locus,show_all_snps = T, drop_hilit = 0.15*7.048437)


#get the colocalizing eqtl
load("./results/Rdata/local_eqtl.Rdata")#eqtl

x = local_eqtl[which(local_eqtl$Gene.Name %in% genes_locus$Name),"lodcolumn"]
x = paste0(x,"_3")

load("~/Desktop/merge_analysis/merge_top_local_eqtl.Rdata") #merge frame

merge_top_maar = merge_top[which(names(merge_top) %in% x)]

coloc_snps = c()
for(i in 1:length(merge_top_maar)){
  if(any(top_maar3$snp_id %in% merge_top_maar[[i]]$snp_id)){
    print(names(merge_top_maar)[i])
    coloc_snps = append(coloc_snps,top_maar3[which(top_maar3$snp_id %in% merge_top_maar[[i]]$snp_id),"snp_id"])
  }
}


coloc_snps = unique(coloc_snps) #colocalizing snps

#create_variant_query_func("./data/CCdb/cc_variants.sqlite", filter="snp_id == " ? ")

#l = out_snps_maar_3$lod[which(rownames(out_snps_maar_3$lod) %in% coloc_snps),]
#s = out_snps_maar_3$snpinfo[which(out_snps_maar_3$snpinfo$snp_id %in% coloc_snps),]




#plot_genes
g = genes_locus[which(genes_locus$start >= start & genes_locus$stop <= end),]
plot_genes(genes = g)

####
####
#Figure3
####
####3A

plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD",ylim=c(0,25),col="red")
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "ML",col="blue",add=T)
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI",col="green",add=T)
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax",col="black",add=T)
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar",col="violet",add=T)
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar",col="orange",add=T)
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar",col="yellow",add=T)
plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity",col="purple",add=T)


patchwork = wrap_plots(x1,x2,x3,x4,x5,x6,x7,x8)
patchwork * theme(axis.text.x = element_blank(),axis.ticks = element_blank(), axis.title = element_blank())

qtl = read.csv("./results/flat/qtl_loc",stringsAsFactors = F)
qtl$pheno_name = c("ML","Ma.Ar","Tt.Ar","TMD","Ct.Por","pMOI","Imax","Ct.Ar/Tt.Ar","MAT_VOL1","ML","Ma.Ar","Ma.Ar","Tt.Ar","Ct.Ar/Tt.Ar","BMD","Disp. @ frax","Disp. @ max load","Tot.Work","WPY","TMD","Max Load","Frax. Load","Ct.Ar","pMOI","Imax","Imin","Tb.Sp","Tb.N","Ct.Th")
qtl = qtl[which(qtl$chr=="1"),]
qtl$chr = as.numeric(qtl$chr)



ggplot() + geom_point(data = qtl, aes(x=pos,y=lod))+ labs(x = "Chromosome 1", y = "LOD") +
  geom_text_repel(data=qtl,aes(x=pos,y=lod,label=pheno_name),size = 3.5)+
  theme(legend.position = "none",panel.grid = element_blank())



###
#plot as lines only?
dat = as.data.frame(DO_qtl_scan_normal[,c("uCT_Ct.TMD","ML","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar","uCT_Tt.Ar","uCT_Ma.Ar","uCT_Ct.porosity")])
map = as.data.frame(cross_basic$pmap$`1`)
dat = as.data.frame(dat[which(rownames(dat) %in% rownames(map)),])
colnames(dat) = c("TMD","ML","pMOI","Imax","Ct.Ar/Tt.Ar","Tt.Ar","Ma.Ar","Por")                  
dat$map = map[1:nrow(dat),1]

x1 = ggplot(data=dat[5700:6200,], aes(x=map, y=TMD))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 8.07,col="red")
x2 = ggplot(data=dat[5700:6200,], aes(x=map, y=ML))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.76,col="red")
x3 = ggplot(data=dat[5700:6200,], aes(x=map, y=pMOI))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.80,col="red")
x4 = ggplot(data=dat[5700:6200,], aes(x=map, y=Imax))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.86,col="red")
x5 = ggplot(data=dat[5700:6200,], aes(x=map, y=`Ct.Ar/Tt.Ar`))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.59,col="red")
x6 = ggplot(data=dat[5700:6200,], aes(x=map, y=Tt.Ar))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.72,col="red")
x7 = ggplot(data=dat[5700:6200,], aes(x=map, y=Ma.Ar))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.63,col="red")
x8 = ggplot(data=dat[5700:6200,], aes(x=map, y=Por))+geom_line()+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+geom_hline(yintercept = 7.77,col="red")+xlab("Chromosome 1 Position")

x1/x2/x3/x4/x5/x6/x7/x8

##3B
#POMP data

#3C 
#like 3A but after conditioning on ML index variant

#condition on top snp rs50769082
snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.4623),
                      sdp=128,
                      snp=c("rs50769082"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs50769082[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs50769082[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus1_scan_cond = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)

dat = as.data.frame(locus1_scan_cond[,c("uCT_Ct.TMD","ML","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar","uCT_Tt.Ar","uCT_Ma.Ar","uCT_Ct.porosity")])
map = as.data.frame(cross_basic$pmap$`1`)
dat = as.data.frame(dat[which(rownames(dat) %in% rownames(map)),])
colnames(dat) = c("TMD","ML","pMOI","Imax","Ct.Ar/Tt.Ar","Tt.Ar","Ma.Ar","Por")                  
dat$map = map[1:nrow(dat),1]

x1 = ggplot(data=dat[5700:6200,], aes(x=map, y=TMD))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 8.07,col="red")
x2 = ggplot(data=dat[5700:6200,], aes(x=map, y=ML))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.76,col="red")
x3 = ggplot(data=dat[5700:6200,], aes(x=map, y=pMOI))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.80,col="red")
x4 = ggplot(data=dat[5700:6200,], aes(x=map, y=Imax))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.86,col="red")
x5 = ggplot(data=dat[5700:6200,], aes(x=map, y=`Ct.Ar/Tt.Ar`))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.59,col="red")
x6 = ggplot(data=dat[5700:6200,], aes(x=map, y=Tt.Ar))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.72,col="red")
x7 = ggplot(data=dat[5700:6200,], aes(x=map, y=Ma.Ar))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.63,col="red")
x8 = ggplot(data=dat[5700:6200,], aes(x=map, y=Por))+geom_line()+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+geom_hline(yintercept = 7.77,col="red")+xlab("Chromosome 1 Position")

x1/x2/x3/x4/x5/x6/x7/x8


