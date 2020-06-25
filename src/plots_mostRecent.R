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
library(tidyr)
library(reshape2)
library(emmeans)
library(ggsignif)
library(RACER)
library(PhenStat)
library(Seurat)
library(igraph)
library(bnlearn)
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
##

#create a covar object for eqtl from covariates in cross file
#must be numeric
covar_eqtl = as.matrix(cross_eqtl$covar)
covar_eqtl[,"sex"] = (covar_eqtl[,"sex"] == "M")*1 #convert sex to 1's and 0's
covar_eqtl[,1] = as.factor(covar_eqtl[,1]) #sac date to factors
covar_eqtl[,6] = as.factor(covar_eqtl[,6]) #generation to factors

covar_eqtl = apply(covar_eqtl,2,as.numeric)
rownames(covar_eqtl) = rownames(cross_eqtl$covar)




##########
#mat is never normal, heavily right skewed, even after adding 1 to remove the exclusion of zero counts
#for pyd, herit is higher if you add 1?
norm_pheno = as.data.frame(cross_basic$pheno)

norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1+1
norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2+1
norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1

norm_pheno$bending_work_post_yield = norm_pheno$bending_work_post_yield + 1
norm_pheno$bending_PYD = norm_pheno$bending_PYD + 1

norm_pheno = as.data.frame(log10(norm_pheno[,c(6:14,16,17,21,23:33,34,35,37:41,43:49,51,52,54:58,61:70,72,74,76)]))

pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,36,42,50,53,59,60)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}



new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

##########


#1B)
## allele freq per geno
#calculate allele freq of each marker, get chromosome from find_markerpos, plot contribution of allele per chrom
chr = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
af_mar <- calc_geno_freq(apr, "marker", omit_x=FALSE)
af_mar = as.data.frame(af_mar)
af_mar$chrom = NA
af_mar$markerName = rownames(af_mar)


af_mar$chrom = as.character(apply(af_mar,1,function(x) find_markerpos(cross_basic$pmap,(x[10]))[1][[1]]))


mean_per_chrom = aggregate(af_mar[,c(1:8)],by = list(af_mar$chrom), FUN = mean)

mean_per_chrom = mean_per_chrom[match(chr,mean_per_chrom$Group.1),]
rownames(mean_per_chrom) =mean_per_chrom$Group.1
mean_per_chrom = mean_per_chrom[,-1]

cairo_pdf(file="~/Desktop/figs/S1A.pdf", width = 10, height = 7)
barplot(as.matrix(t(mean_per_chrom)),col = CCcolors, ylab = "Allele frequency", xlab="Chr", main = "Global allele frequency per chromosome")
#legend("topleft", fill=CCcolors, legend=c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ","NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ"),cex = 2, box.lwd = 1.5)
dev.off()





# dat = mean_per_chrom %>% gather()
# dat$chrom = rep(rownames(mean_per_chrom), nrow(dat)/nrow(mean_per_chrom))
# 
# dat$chrom = as.factor(dat$chrom)
# levels(dat$chrom) = c(1:20,"X")
# ggplot(dat, aes(x=chrom, y= value, fill=key)) + geom_bar(position=position_fill(reverse=T ),stat = "identity")
# 
# p+ scale_fill_manual(values=as.vector(CCcolors))

# 
# par(mar=c(4.1, 4.1, 0.6, 0.6))
# 
# par(mfrow=c(4,2))
# for(i in 1:8){ hist(af_ind[,i], main=NULL, breaks=30,xlab = "geno freq by ind.",cex.lab=1.25, cex.axis=1.25)
#   abline(v=mean(af_ind[,i]),col="red",lwd=3)
# }



#1C)
#boxplot BV/TV
#make dataframe with bv/tv 
df_bvtv = as.data.frame(cross_basic$pheno[,"uCT_BV.TV"])
df_bvtv$pheno = "BV/TV"
df_bvtv$ind = rownames(df_bvtv)
colnames(df_bvtv)[1] = "val"


p3 <- ggplot(df_bvtv, aes(x=reorder(ind,val), y=val,width=1)) + 
  geom_bar(stat="identity",color=CCcolors[4]) + theme_bw() + theme(axis.text.x = element_blank()) + xlab("") + ylab("Bone volume fraction (BV/TV, %)") + xlab("DO mouse index") + scale_y_continuous(expand=c(0,0)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_text(size=20), axis.text.y = element_text(size = 15))


df_str = as.data.frame(cross_basic$pheno[,"bending_max_load"])
df_str$pheno = "max_load"
df_str$ind = rownames(df_str)
colnames(df_str)[1] = "val"

p4 <- ggplot(df_str, aes(x=reorder(ind,val), y=val,width=1)) + 
  geom_bar(stat="identity",color=CCcolors[6]) +theme_bw() + theme(axis.text.x = element_blank()) + xlab("DO Mouse") + ylab("Max Load (N)") + xlab("DO mouse index") + scale_y_continuous(expand=c(0,0)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_text(size=20), axis.text.y = element_text(size = 15))





patch = p3 / p4

pdf(file="~/Desktop/figs/1C.pdf", width = 10, height = 7)
patch
dev.off()

cairo_pdf(file="~/Desktop/figs/1Ca.pdf", width = 10, height = 7)
p3
dev.off()

cairo_pdf(file="~/Desktop/figs/1Cb.pdf", width = 10, height = 7)
p4
dev.off()
# 
# p
#1D) heritability


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
hsq_table = h_df[-which(h_df$pheno %in% c("glucose", "soleus_weight","adiposity","RFP","GFP","BFP","FFP","gastroc_weight","MAT_VOL1_nonzero","MAT_VOL2_nonzero","MAT_VOL3_nonzero","MAT_VOL4_nonzero","body_length")),]

hsq_table[grep("bending",x = hsq_table$pheno),"group"] = 1
hsq_table[grep("uCT",x = hsq_table$pheno),"group"] = 2
hsq_table[grep("histo",x = hsq_table$pheno),"group"] = 3
hsq_table[grep("MAT",x = hsq_table$pheno),"group"] = 4

hsq_table[grep("ML",x = hsq_table$pheno,fixed = T),"group"] = 5
hsq_table[grep("AP",x = hsq_table$pheno,fixed=T),"group"] = 5
hsq_table[grep("FL",x = hsq_table$pheno,fixed=T),"group"] = 5
hsq_table[grep("body_length",x = hsq_table$pheno,fixed=T),"group"] = 5

hsq_table = hsq_table[order(hsq_table$h),]


pheno = c(hsq_table[which(hsq_table$group == 1),"pheno"])
pheno = append(pheno,c(hsq_table[which(hsq_table$group == 2),"pheno"]))
pheno = append(pheno,c(hsq_table[which(hsq_table$group == 3),"pheno"]))
pheno = append(pheno,c(hsq_table[which(hsq_table$group == 4),"pheno"]))
pheno = append(pheno,c(hsq_table[which(hsq_table$group == 5),"pheno"]))







hsq_table$pheno <- factor(hsq_table$pheno, levels = hsq_table$pheno)



p5<-ggplot(hsq_table, aes(x=pheno, y=h, fill=as.factor(group))) + 
  geom_bar(stat="identity",show.legend = F) + coord_flip() + scale_x_discrete(limits = c(pheno), labels=c(
    expression(D[yield]),"PYD", expression(W[yield]), "k", expression(k[yield]),
    expression(F[yield]),expression(D[Fmax]), expression(W[py]), "W", expression(D[fx]),
    expression(F[max]),expression(F[fx]), "Ct.Por","SMI","Tb.Th",
    "BS/BV","Conn.D.", "BV/TV","BMD","Tb.Sp",
    "Tb.N", expression(I[max]), "TMD",expression(I[min]),"pMOI",
    "Tt.Ar","Ct.Th", "Ct.Ar","Ma.Ar","Ct.Ar/Tt.Ar",
    "histo_OS/BS", "histo_OV/BV", "histo_ObS/OS", "histo_NOb/Obpm", "histo_NOb/Opm",
    "histo_BS/BV", "histo_Tb.Th", "histo_NOc/Bpm", "histo_OTh", "histo_ObS/BS",
    "histo_OcS/BS", "histo_Tb.N", "histo_NOb/Bpm", "histo_BV/TV", "histo_Tb.Sp",
    "histo_NOc/TAR", "histo_NOb/TAR","histo_OV/TV", "MAT_V1", "MAT_V3",
    "MAT_V2", "MAT_V4","AP","FL","ML"
    )) + 
  
  xlab("Phenotype") + ylab("Heritability") + scale_color_brewer(palette = "Dark2") + theme_bw() +scale_y_continuous(limits = c(0,1),expand=c(0,0)) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_text(size=20), axis.text = element_text(size = 15))

cairo_pdf(file="~/Desktop/figs/S1B.pdf", width = 12, height = 14)
p5
dev.off()

s#wrap_plots(p1,p,p2,p3)
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


p2 = wrap_elements(panel = ~barplot(as.matrix(t(mean_per_chrom)),col = CCcolors, ylab = "Allele Frequency", xlab="Chr", main = "Global Allele Freq. Per Chromosome"), clip=FALSE)

patch = wrap_elements(patch)

layout = "
AABB
CCDD
"

p1+p1+patch+p5+plot_layout(design = layout)+plot_annotation(tag_levels = "A")



#
#
#
#

####################################################################################################################################################################################
####################################################################################################################################################################################
####################################################################################################################################################################################

##Figure 2
#2A - mirrorplot
#Sertad4 - adipose_subcutaneous

morris_lead_snps = read.csv("./data/Morrisetal2018.NatGen.SumStats/Morris_eBMD_conditionally_ind_snps.csv", header=T, stringsAsFactors = F)

ebmd = read.delim("./data/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("./data/adipose_subcutaneous_SERTAD4_GTEX_v7", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]


#gtex$a1 =  sapply(strsplit(gtex$V2, "_"), "[[", 3)
#gtex$a2 =  sapply(strsplit(gtex$V2, "_"), "[[", 4)

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]

#remove duplicated by matching alleles in ebmd and gtex
ebmd[which(duplicated(ebmd$RSID)), "RSID"]
 #317,636,907, 1099
ebmd = ebmd[-c(317,636,907,1099),]







for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"pos"]


bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 15, p_col = 13)
bmd_racer = bmd_racer[,-11]
#


gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs7516171"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs7516171"))

cairo_pdf(file="~/Desktop/figs/2A_sertad.pdf", width = 9, height = 7)
par(ps=12)
mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 1, plotby = "snp",snp_plot = "rs7516171",build = "hg19", name1 = "eBMD", name2 = "SERTAD4 - Adipose Subcutaneous")
mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 1, plotby = "coord",start_plot = 210100000,end_plot = 210850000,build = "hg19", name1 = "eBMD", name2 = "SERTAD4 - Adipose Subcutaneous")

dev.off()
####


####################################################################################################################################################################################
#Glt8d2 pituitary

morris_lead_snps = read.csv("./data/Morrisetal2018.NatGen.SumStats/Morris_eBMD_conditionally_ind_snps.csv", header=T, stringsAsFactors = F)

ebmd = read.delim("./data/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("./data/pituitary_GLT8D2_GTEX_V7", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]


#gtex$a1 =  sapply(strsplit(gtex$V2, "_"), "[[", 3)
#gtex$a2 =  sapply(strsplit(gtex$V2, "_"), "[[", 4)

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]

#remove duplicated by matching alleles in ebmd and gtex
ebmd[which(duplicated(ebmd$RSID)), "RSID"]
#71 642 942 981 988
ebmd = ebmd[-c(71, 642, 942, 981, 988),]







for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"pos"]


bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 15, p_col = 13)
bmd_racer = bmd_racer[,-11]
#


gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs2722176"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs2722176"))

cairo_pdf(file="~/Desktop/figs/2A_glt8d2.pdf", width = 9, height = 7)
par(ps=12)
mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 12, plotby = "coord", start_plot = 103900000, end_plot = 104550000, build = "hg19", name1 = "eBMD", name2 = "GLT8D2 - Pituitary")
dev.off()
####

#2B
net = readRDS("./results/Rdata/networks/wgcna_m_5.RDS")
rnames = rownames(net$MEs)
royalblue = as.data.frame(net$MEs[,"ME20"])
royalblue$DO = rnames
royalblue = royalblue[order(as.numeric(royalblue$DO)),]

tbsp = pheno_combined[,"uCT_Tb.Sp"]

n=c()
for(i in names(tbsp)){
  x=strsplit(i, split = "[.]")[[1]][1]
  n=append(n,x)
}
names(tbsp) = n
tbsp = tbsp[which(names(tbsp) %in% royalblue$DO)]

tbsp = tbsp[order(as.numeric(names(tbsp)))]

cor(royalblue$`net$MEs[, "ME20"]`, tbsp, method = "s", use = "p")
cor.test(royalblue$`net$MEs[, "ME20"]`, tbsp,method = "s", use="p")


d = cbind(tbsp,royalblue)
colnames(d)[2] = "ME20"
#rho 0.2732899
#pval 0.007058
pdf(file="~/Desktop/figs/2B1.pdf", width = 9, height = 7)
ggplot(d, aes(y=ME20, x=tbsp)) + geom_point() + geom_smooth(method=lm) + ylab("Royalblue_m module eigengenes") + xlab("Trabecular Separation ()") + 
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_text(size=18), axis.text = element_text(size = 15))

dev.off()








p = pheno_combined[,"uCT_Tb.N"]

n=c()
for(i in names(p)){
  x=strsplit(i, split = "[.]")[[1]][1]
  n=append(n,x)
}
names(p) = n
p = p[which(names(p) %in% royalblue$DO)]

p = p[order(as.numeric(names(p)))]

cor(royalblue$`net$MEs[, "ME20"]`, p, method = "s", use = "p") #-0.263
cor.test(royalblue$`net$MEs[, "ME20"]`, p,method = "s", use="p")#p=0.009517


d = cbind(p,royalblue)
colnames(d)[2] = "ME20"


pdf(file="~/Desktop/figs/2B2.pdf", width = 9, height = 7)
ggplot(d, aes(y=ME20, x=p)) + geom_point() + geom_smooth(method=lm) + ylab("Roaylblue module eigengenes") + xlab(expression(paste("Trabecular Number "))) +
  theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_text(size=18), axis.text = element_text(size = 15))

dev.off()




# layout = "
# ABC
# DEF
# "
# 
# par(ps=12)
# cairo_pdf(file="~/Desktop/figs/2B.pdf", width = 15, height = 12)
# p1+p2+p3+p4+p5+p6+plot_layout(design = layout)
# dev.off()
# 

####################################################################################################################################################################################

#2C
bn2igraph <- function(g.bn){
  g <- igraph.from.graphNEL(as.graphNEL(g.bn))
}




load("./results/Rdata/networks/bn_m_5/hybrid_royalblue_m_5.RData")

assign(x = "royalblue_bn" ,bn)
#convert and plot a particular BN
x = bn2igraph(royalblue_bn)


subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,2,nodes = "Glt8d2"))))
plot(subgraph,vertex.label.cex=0.65,edge.width=2, vertex.size=10, margin=-0.1, vertex.label.dist=0.2, vertex.label.degree=-pi)


subgraph <- induced.subgraph(x, names(unlist(neighborhood(x,2,nodes = "Sertad4"))))
plot(subgraph,vertex.label.cex=0.65,edge.width=2, vertex.size=10, margin=-0.1, vertex.label.dist=0.2, vertex.label.degree=-pi)


  ####################################################################################################################################################################################
#2D
#FROM B6_OBs_RNA_SEQ
####################################################################################################################################################################################

#2E
load("./results/Rdata/seurat_ob.Rdata") #loads as "ob

cairo_pdf(file="~/Desktop/figs/qsox1_vln.pdf", width = 10, height = 7)
VlnPlot(ob,features = "Qsox1")
dev.off()

cairo_pdf(file="~/Desktop/figs/glt8d2_vln.pdf", width = 10, height = 7)
VlnPlot(ob,features = "Glt8d2")
dev.off()

cairo_pdf(file="~/Desktop/figs/sertad4_vln.pdf", width = 10, height = 7)
VlnPlot(ob,features = "Sertad4")
dev.off()

cairo_pdf(file="~/Desktop/figs/2E_umap.pdf", width = 10, height = 7)
DimPlot(ob,
        reduction = "umap",
        label = TRUE,
        label.size = 6,
        pt.size = 1.5)
dev.off()

cairo_pdf(file="~/Desktop/figs/2E_sertad4.pdf", width = 10, height = 7)
FeaturePlot(ob, features = c("Sertad4"),pt.size = 1.5, sort.cell = T)
dev.off()



cairo_pdf(file="~/Desktop/figs/2E_Glt8d2.pdf", width = 10, height = 7)
FeaturePlot(ob, features = c("Glt8d2"),pt.size = 1.5, sort.cell = T)
dev.off()

####################################################################################################################################################################################
#2F
library(IMPCdata)
t1<-getIMPCDataset('WTSI','MGP_001','IMPC_DXA_001','IMPC_DXA_004_001','MGI:4364018')
library(PhenStat)
t1<-subset(t1,t1$Genotype!="colony_id")
plot(as.vector(as.character(t1$Weight)),as.numeric(as.character(t1$Value)))
t1$Value<-as.numeric(as.character(t1$Value))
t1$Weight<-as.numeric(as.character(t1$Weight))
o1<-unique(t1$Genotype)
d2<-PhenList(t1, testGenotype=as.character(unique(t1$Genotype)[which(o1!='+/+')]),refGenotype='+/+',outputMessages = F)
d3<-testDataset(d2,depVariable='Value',method='MM',equation='withWeight', outputMessages = T,transformValues = T)
#d3@analysisResults$model.output.summary.

cairo_pdf(file="~/Desktop/figs/2F.pdf", width = 10, height = 7)
ggplot(t1,aes(x=Sex, y=Value, fill=factor(Genotype))) +
  geom_boxplot() +
  labs(fill = "Genotype") + 
  geom_point(position=position_jitterdodge(),alpha=0.3) +
  theme_bw(base_size = 16)
dev.off()
# ##2D
# rasd1_impc = read.csv("~/Desktop/test.txt")
# #r_m = rasd1_impc[which(rasd1_impc$Sex == "male"),]
# #r_f = rasd1_impc[which(rasd1_impc$Sex == "female"),]
# r = rasd1_impc
# 
# r$con = paste0(r$Genotype, r$Sex)
# r[which(r$con == "+/+female"),"con"] = "Female WT"
# r[which(r$con == "+/+male"),"con"] = "Male WT"
# r[which(r$con == "BL3486female"),"con"] = "KO Female"
# r[which(r$con == "BL3486male"),"con"] = "KO Male"
# 
# r$con = factor(r$con, levels = c("Female WT","KO Female","Male WT", "KO Male"))
# r$col = 1
# r$col[which(r$con == "Female WT" | r$con == "Male WT")] = 0
# x = PhenList(rasd1_impc, testGenotype = "EPD0098_5_B05",refGenotype = "+/+")
# 
# t = testDataset(x, depVariable = "Value", equation = "withWeight")
# 
# par(ps=12)
# cairo_pdf(file="~/Desktop/figs/2D.pdf", width = 10, height = 7)
# ggplot(r, aes(con,Value, fill=as.factor(col) )) + geom_boxplot(width=0.1) + xlab(NULL) + ylab("BMD") + theme_bw() + theme(panel.grid = element_blank(), panel.border = element_blank(), axis.title = element_text(size=18), axis.text = element_text(size = 15), legend.position = "none")+
#   scale_fill_brewer(palette = "Dark2")
# dev.off()
# 
# #resid = as.data.frame(t@analysisResults$model.output$residuals)
# 
# #resid$genotype = rasd1_impc[,"Genotype"]
# #resid$sex = rasd1_impc[,"Sex"]
# 
# ##2C
# load("./results/Rdata/seurat_ob.Rdata") #loads as "ob
# 
# #plot UMAP
# ob=JackStraw(ob)
# ob <- RunUMAP(ob, dims = 1:13)
# 
# cairo_pdf(file="~/Desktop/figs/2C_umap.pdf", width = 10, height = 7)
# DimPlot(ob,
#         reduction = "umap",
#         label = TRUE,
#         label.size = 6,
#         plot.title = "UMAP")
# dev.off()
# 
# cairo_pdf(file="~/Desktop/figs/2C_rasd1.pdf", width = 10, height = 7)
# FeaturePlot(ob, features = c("Rasd1"),pt.size = 1.5, sort.cell = T)
# dev.off()
# 
# 
# 
# ob.markers <- FindAllMarkers(ob, only.pos = TRUE)
# 
# 
# #First, Rasd1 markers. (Cluster 10)
# 
# ob.markers[which(ob.markers$gene =="Rasd1"),] # 10
# 
# #output the top 30 cluster 10 genes
# cluster10.markers <- FindMarkers(ob, ident.1 = 10, min.pct = 0.25,)
# 
# head(cluster10.markers, 30)
# 
# 
# cairo_pdf(file="~/Desktop/figs/2C_sost.pdf", width = 10, height = 7)
# FeaturePlot(ob, features = c("Sost"),pt.size = 1.5, sort.cell = T)
# dev.off()
# 
# cairo_pdf(file="~/Desktop/figs/2C_phex1.pdf", width = 10, height = 7)
# FeaturePlot(ob, features = c("Phex"),pt.size = 1.5, sort.cell = T)
# dev.off()

##2E
#from B6>OBs_RNA_Seq.R
# load("./results/Rdata/rasd1_calvarial.Rdata") #called temp4
# 
# genes = "Rasd1"
# 
# pdf(file="~/Desktop/figs/2E.pdf", width = 10, height = 7)
# ggplot(temp4, aes(x=Day, y=mean,ymin=mean-temp4$sem, ymax=mean+temp4$sem, colour=newGene,Group=newGene)) + 
#   geom_errorbar(aes(ymin=mean-temp4$sem, ymax=mean+temp4$sem), width=0.5, size=1.25) +
#   geom_line(aes(group=newGene), size=0.5) +
#   geom_point() + scale_color_discrete(name="Gene",
#                                       breaks=as.character(unique(temp4$Gene)),
#                                       labels=genes,l=60) +
#   scale_y_continuous() + facet_wrap(~newGene, scales='free_y') 
# dev.off()
# 







##Figure 3
#3A
#scatterplot with traits, across all chroms, color traits differently, LOD score, labels
#load in significant qtl
qtl = read.csv("./results/flat/qtl_loc",stringsAsFactors = F)
qtl$pheno_name = c("ML","Ma.Ar","Tt.Ar","TMD","Ct.Por","pMOI","Imax","Ct.Ar/Tt.Ar","ML","Ma.Ar","Ma.Ar","Tt.Ar","Ct.Ar/Tt.Ar","BMD","Dfx","DFmax","Wtotal","WPY","TMD","Fmax","Ffx","Ct.Ar","pMOI","Imax","Imin","Tb.Sp","Tb.N","Ct.Th")
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


p2a=ggplot() + 
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

cairo_pdf(file="~/Desktop/figs/S2.pdf", width = 10, height = 7)
p2a
dev.off()







###4A
#chr3 Ma.Ar QTL Trace

#load("./results/Rdata/DO_qtl_scan_norm.Rdata")

#maar3_blup = scan1blup(apr[,3], pheno_combined[,"uCT_Ma.Ar"], k_loco[["3"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
#save(maar3_blup,file = "./results/Rdata/maar_chr3_blup.Rdata")

load("./results/Rdata/maar_chr3_blup.Rdata")
cairo_pdf(file="~/Desktop/figs/4A.pdf", width = 10, height = 7)
plot_coefCC(maar3_blup[1000:3300,], cross_basic$pmap,scan1_output = subset(DO_qtl_scan_normal, lodcolumn="uCT_Ma.Ar"), main = "Ma.Ar - Chr. 3", legend_ncol=1,top_panel_prop = 0.6)
dev.off()




##4B

query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")



start = 66.56421 - 0.250
end = 69.96983 + 0.250
chr = 3
out_snps_maar_3 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ma.Ar"], k_loco[["3"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                           query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=T)


variants_locus = query_variants(chr, start, end)
genes_locus <- query_genes(chr, start, end)

#genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]

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

#pdf(file="~/Desktop/figs/3C.pdf", width = 10, height = 7)
#plot_snpasso(out_snps_maar_3$lod, out_snps_maar_3$snpinfo, genes = genes_locus,show_all_snps = T, drop_hilit = 0.15*7.048437)
#dev.off()

#get the colocalizing eqtl
load("./results/Rdata/local_eqtl.Rdata")#eqtl

x = local_eqtl[which(local_eqtl$Gene.Name %in% genes_locus$Name),"lodcolumn"]
x = paste0(x,"_3")

load("./results/Rdata/merge_top_local_eqtl.Rdata") #merge frame

merge_top_maar = merge_top[which(names(merge_top) %in% x)]

coloc_snps = c()
for(i in 1:length(merge_top_maar)){
  if(any(top_maar3$snp_id %in% merge_top_maar[[i]]$snp_id)){
    print(names(merge_top_maar)[i])
    coloc_snps = append(coloc_snps,top_maar3[which(top_maar3$snp_id %in% merge_top_maar[[i]]$snp_id),"snp_id"])
  }
}

#coloc genes: MSTRG.15696, MSTRG.15702, MSTRG.15703, MSTRG.15704
load("./results/Rdata/chr3_coloc_genes.Rdata") #called "chr3_coloc_genes", got from supercomputing cluster. files too big to keep on laptop
#plot
MSTRG.15696 = scan1blup(apr[,3], cross_eqtl$pheno[,"MSTRG.15696"], kinship = k_loco[[3]], addcovar = covar_eqtl[,c(2,11:58)])
MSTRG.15702 = scan1blup(apr[,3], cross_eqtl$pheno[,"MSTRG.15702"], kinship = k_loco[[3]], addcovar = covar_eqtl[,c(2,11:58)])
MSTRG.15703 = scan1blup(apr[,3], cross_eqtl$pheno[,"MSTRG.15703"], kinship = k_loco[[3]], addcovar = covar_eqtl[,c(2,11:58)])
MSTRG.15704 = scan1blup(apr[,3], cross_eqtl$pheno[,"MSTRG.15704"], kinship = k_loco[[3]], addcovar = covar_eqtl[,c(2,11:58)])


cairo_pdf(file="~/Desktop/figs/4B1.pdf", width = 10, height = 7)
plot_coefCC(MSTRG.15696[1000:3300,], cross_basic$pmap,scan1_output = subset(chr3_coloc_genes, lodcolumn="MSTRG.15696"), main = "Mfsd1 - Chr. 3", legend_ncol=1,top_panel_prop = 0.6)
dev.off()

cairo_pdf(file="~/Desktop/figs/4B2.pdf", width = 10, height = 7)
plot_coefCC(MSTRG.15702[1000:3300,], cross_basic$pmap,scan1_output = subset(chr3_coloc_genes, lodcolumn="MSTRG.15702"), main = "Il12a - Chr. 3", legend_ncol=1,top_panel_prop = 0.6)
dev.off()

cairo_pdf(file="~/Desktop/figs/4B3.pdf", width = 10, height = 7)
plot_coefCC(MSTRG.15703[1000:3300,], cross_basic$pmap,scan1_output = subset(chr3_coloc_genes, lodcolumn="MSTRG.15703"), main = "Gm17641 - Chr. 3", legend_ncol=1,top_panel_prop = 0.6)
dev.off()

cairo_pdf(file="~/Desktop/figs/4B4.pdf", width = 10, height = 7)
plot_coefCC(MSTRG.15704[1000:3300,], cross_basic$pmap,scan1_output = subset(chr3_coloc_genes, lodcolumn="MSTRG.15704"), main = "1110032F04Rik - Chr. 3", legend_ncol=1,top_panel_prop = 0.6)
dev.off()
# 




























####
####

####
####5A

# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD",ylim=c(0,25),col="red")
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "ML",col="blue",add=T)
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI",col="green",add=T)
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax",col="black",add=T)
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar",col="violet",add=T)
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar",col="orange",add=T)
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar",col="yellow",add=T)
# plot_scan1(DO_qtl_scan_normal[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity",col="purple",add=T)
# 
# 
# patchwork = wrap_plots(x1,x2,x3,x4,x5,x6,x7,x8)
# patchwork * theme(axis.text.x = element_blank(),axis.ticks = element_blank(), axis.title = element_blank())
# 
# qtl = read.csv("./results/flat/qtl_loc",stringsAsFactors = F)
# qtl$pheno_name = c("ML","Ma.Ar","Tt.Ar","TMD","Ct.Por","pMOI","Imax","Ct.Ar/Tt.Ar","MAT_VOL1","ML","Ma.Ar","Ma.Ar","Tt.Ar","Ct.Ar/Tt.Ar","BMD","Disp. @ frax","Disp. @ max load","Tot.Work","WPY","TMD","Max Load","Frax. Load","Ct.Ar","pMOI","Imax","Imin","Tb.Sp","Tb.N","Ct.Th")
# qtl = qtl[which(qtl$chr=="1"),]
# qtl$chr = as.numeric(qtl$chr)
# 
# 
# 
# ggplot() + geom_point(data = qtl, aes(x=pos,y=lod))+ labs(x = "Chromosome 1", y = "LOD") +
#   geom_text_repel(data=qtl,aes(x=pos,y=lod,label=pheno_name),size = 3.5)+
#   theme(legend.position = "none",panel.grid = element_blank())
# 


###
#plot as lines only
dat = as.data.frame(DO_qtl_scan_normal[,c("uCT_Ct.TMD","ML","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar","uCT_Tt.Ar","uCT_Ma.Ar","uCT_Ct.porosity")])
map = as.data.frame(cross_basic$pmap$`1`)
dat = as.data.frame(dat[which(rownames(dat) %in% rownames(map)),])
colnames(dat) = c("TMD","ML","pMOI","Imax","Ct.Ar/Tt.Ar","Tt.Ar","Ma.Ar","Por")                  
dat$map = map[1:nrow(dat),1]

# x1 = ggplot(data=dat[5700:6200,], aes(x=map, y=TMD))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 8.07,col="red")
# x2 = ggplot(data=dat[5700:6200,], aes(x=map, y=ML))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.76,col="red")
# x3 = ggplot(data=dat[5700:6200,], aes(x=map, y=pMOI))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.80,col="red")
# x4 = ggplot(data=dat[5700:6200,], aes(x=map, y=Imax))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.86,col="red")
# x5 = ggplot(data=dat[5700:6200,], aes(x=map, y=`Ct.Ar/Tt.Ar`))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.59,col="red")
# x6 = ggplot(data=dat[5700:6200,], aes(x=map, y=Tt.Ar))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.72,col="red")
# x7 = ggplot(data=dat[5700:6200,], aes(x=map, y=Ma.Ar))+geom_line()+theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.x = element_blank())+geom_hline(yintercept = 7.63,col="red")
# x8 = ggplot(data=dat[5700:6200,], aes(x=map, y=Por))+geom_line()+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+geom_hline(yintercept = 7.77,col="red")+xlab("Chromosome 1 Position")
# 
# cairo_pdf(file="~/Desktop/figs/5A.pdf", width = 10, height = 7)
# x1/x2/x3/x4/x5/x6/x7/x8
# dev.off()

TMD_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.TMD"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
ML_blup = scan1blup(apr[,1], pheno_combined[,"ML"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
pMOI_blup = scan1blup(apr[,1], pheno_combined[,"uCT_pMOI"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
Imax_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Imax"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
ctarttar_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.Ar.Tt.Ar"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
ttar_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Tt.Ar"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
maar_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Ma.Ar"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
por_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.porosity"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)

#save(TMD_blup, ML_blup, pMOI_blup, Imax_blup, ctarttar_blup, ttar_blup, maar_blup, por_blup, file = "~/Desktop/5Ablups.Rdata")
cairo_pdf(file="~/Desktop/figs/5Ablups1.pdf", width = 10, height = 7)
plot_coefCC(TMD_blup[5700:6200,], cross_basic$pmap, main = "TMD - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()

cairo_pdf(file="~/Desktop/figs/5Ablups2.pdf", width = 10, height = 7)
plot_coefCC(ML_blup[5700:6200,], cross_basic$pmap, main = "ML - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ablups3.pdf", width = 10, height = 7)
plot_coefCC(pMOI_blup[5700:6200,], cross_basic$pmap, main = "pMOI - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ablups4.pdf", width = 10, height = 7)
plot_coefCC(Imax_blup[5700:6200,], cross_basic$pmap, main = "Imax - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()

cairo_pdf(file="~/Desktop/figs/5Ablups5.pdf", width = 10, height = 7)
plot_coefCC(ctarttar_blup[5700:6200,], cross_basic$pmap, main = "Ct.AR/Tt.Ar - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ablups6.pdf", width = 10, height = 7)
plot_coefCC(ttar_blup[5700:6200,], cross_basic$pmap, main = "Tt.Ar - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ablups7.pdf", width = 10, height = 7)
plot_coefCC(maar_blup[5700:6200,], cross_basic$pmap, main = "Ma.Ar - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ablups8.pdf", width = 10, height = 7)
plot_coefCC(por_blup[5700:6200,], cross_basic$pmap, main = "Ct.Por - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()




cairo_pdf(file="~/Desktop/figs/5Ascan1.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.TMD" )
add_threshold(cross_basic$pmap, thresholdA = 8.07, col="red")
dev.off()

cairo_pdf(file="~/Desktop/figs/5Ascan2.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "ML" )
add_threshold(cross_basic$pmap, thresholdA = 7.76, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ascan3.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_pMOI" )
add_threshold(cross_basic$pmap, thresholdA = 7.80, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ascan4.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Imax" )
add_threshold(cross_basic$pmap, thresholdA = 7.86, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ascan5.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.Ar.Tt.Ar" )
add_threshold(cross_basic$pmap, thresholdA = 7.59, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ascan6.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Tt.Ar" )
add_threshold(cross_basic$pmap, thresholdA = 7.72, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ascan7.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ma.Ar")
add_threshold(cross_basic$pmap, thresholdA = 7.63, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Ascan8.pdf", width = 10, height = 7)
plot_scan1(DO_qtl_scan_normal[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.porosity")
add_threshold(cross_basic$pmap, thresholdA = 7.77, col="red")
dev.off()











# ##4B (now supplementary)
# #POMP data
# #from pomp_mapping_ML.R
# load("./results/Rdata/scan1_pomp.Rdata")
# load("./results/Rdata/coef_ML_pomp_blup.Rdata")
# load("./data/POMP/MM_snps1_pmap.Rdata")
# 
# pdf(file="~/Desktop/figs/4B.pdf", width = 10, height = 7)
# plot_coefCC(coef_ML_pomp_blup, MM_snps1_pmap["1"], scan1_output=subset(scan1_pomp, lodcolumn=2),legend = "topleft")
# dev.off()
# 


#5B
#like 5A but after conditioning on ML index variant

#condition on top snp rs50769082
snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.4623),
                      sdp=128,
                      snp=c("rs50769082"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
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

cairo_pdf(file="~/Desktop/figs/5b.pdf", width = 10, height = 7)
x1/x2/x3/x4/x5/x6/x7/x8
dev.off()



TMD_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.TMD"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
ML_blup_c = scan1blup(apr[,1], pheno_combined[,"ML"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
pMOI_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_pMOI"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
Imax_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_Imax"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
ctarttar_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.Ar.Tt.Ar"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
ttar_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_Tt.Ar"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
maar_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_Ma.Ar"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)
por_blup_c = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.porosity"], k_loco[["1"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)

save(TMD_blup_c, ML_blup_c, pMOI_blup_c, Imax_blup_c, ctarttar_blup_c,ttar_blup_c,maar_blup_c, por_blup_c ,file = "~/Desktop/5Bblups.Rdata")


cairo_pdf(file="~/Desktop/figs/5Bblups1.pdf", width = 10, height = 7)
plot_coefCC(TMD_blup_c[5700:6200,], cross_basic$pmap, main = "TMD - Chr. 1", legend_ncol=1,top_panel_prop = 0.7,)
dev.off()

cairo_pdf(file="~/Desktop/figs/5Bblups2.pdf", width = 10, height = 7)
plot_coefCC(ML_blup_c[5700:6200,], cross_basic$pmap, main = "ML - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bblups3.pdf", width = 10, height = 7)
plot_coefCC(pMOI_blup_c[5700:6200,], cross_basic$pmap,main = "pMOI - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bblups4.pdf", width = 10, height = 7)
plot_coefCC(Imax_blup_c[5700:6200,], cross_basic$pmap, main = "Imax - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bblups5.pdf", width = 10, height = 7)
plot_coefCC(ctarttar_blup_c[5700:6200,], cross_basic$pmap, main = "Ct.AR/Tt.Ar - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bblups6.pdf", width = 10, height = 7)
plot_coefCC(ttar_blup_c[5700:6200,], cross_basic$pmap, main = "Tt.Ar - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bblups7.pdf", width = 10, height = 7)
plot_coefCC(maar_blup_c[5700:6200,], cross_basic$pmap, main = "Ma.Ar - Chr. 1", legend_ncol=1,top_panel_prop = 0.7)
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bblups8.pdf", width = 10, height = 7)
plot_coefCC(por_blup_c[5700:6200,], cross_basic$pmap, legend_ncol=1,top_panel_prop = 0.7)
dev.off()





cairo_pdf(file="~/Desktop/figs/5Bscan1.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.TMD" )
add_threshold(cross_basic$pmap, thresholdA = 8.07, col="red")
dev.off()

cairo_pdf(file="~/Desktop/figs/5Bscan2.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "ML" , ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 7.76, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bscan3.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_pMOI", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.80, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bscan4.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Imax", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.86, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bscan5.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.Ar.Tt.Ar", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.59, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bscan6.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Tt.Ar", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.72, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bscan7.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ma.Ar", ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 7.63, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Bscan8.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.porosity", ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 7.77, col="red")
dev.off()









#condition on top snp rs248974780
snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.0559),
                      sdp=128,
                      snp=c("rs248974780"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}





locus1_scan_cond2 = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)

dat = as.data.frame(locus1_scan_cond2[,c("uCT_Ct.TMD","ML","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar","uCT_Tt.Ar","uCT_Ma.Ar","uCT_Ct.porosity")])
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

cairo_pdf(file="~/Desktop/figs/5c.pdf", width = 10, height = 7)
x1/x2/x3/x4/x5/x6/x7/x8
dev.off()




cairo_pdf(file="~/Desktop/figs/5Cscan1.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.TMD", ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 8.07, col="red")
dev.off()

cairo_pdf(file="~/Desktop/figs/5Cscan2.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "ML" , ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 7.76, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Cscan3.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_pMOI", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.80, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Cscan4.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Imax", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.86, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Cscan5.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.Ar.Tt.Ar", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.59, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Cscan6.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Tt.Ar", ylim=c(0,10) )
add_threshold(cross_basic$pmap, thresholdA = 7.72, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Cscan7.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ma.Ar", ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 7.63, col="red")
dev.off()
cairo_pdf(file="~/Desktop/figs/5Cscan8.pdf", width = 10, height = 7)
plot_scan1(locus1_scan_cond2[5700:6200,], cross_basic$pmap,lodcolumn = "uCT_Ct.porosity", ylim=c(0,10))
add_threshold(cross_basic$pmap, thresholdA = 7.77, col="red")
dev.off()











#NOW 6A
#
load("./results/Rdata/qsox1_ier5.Rdata") # got from supercomputing cluster. main files too big to keep on laptop

#plot
#MSTRG.1301 - Ier5
#MSTRG.1311 - Qsox1

ier5= scan1blup(apr[,1], cross_eqtl$pheno[,"MSTRG.1301"], kinship = k_loco[[1]], addcovar = covar_eqtl[,c(2,11:58)])
qsox1 = scan1blup(apr[,1], cross_eqtl$pheno[,"MSTRG.1311"], kinship = k_loco[[1]], addcovar = covar_eqtl[,c(2,11:58)])

cairo_pdf(file="~/Desktop/figs/5A1.pdf", width = 10, height = 7)
plot_coefCC(ier5[5200:6700,], cross_basic$pmap,scan1_output = subset(x, lodcolumn="Ier5"), main = "Ier5 - Chr. 1", legend_ncol=1,top_panel_prop = 0.6)
dev.off()

cairo_pdf(file="~/Desktop/figs/5A2.pdf", width = 10, height = 7)
plot_coefCC(qsox1[5200:6700,], cross_basic$pmap,scan1_output = subset(x, lodcolumn="Qsox1"), main = "Qsox1 - Chr. 1", legend_ncol=1,top_panel_prop = 0.6)
dev.off()


##
TMD_blup = scan1blup(apr[,1], pheno_combined[,"uCT_Ct.TMD"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)
ML_blup = scan1blup(apr[,1], pheno_combined[,"ML"], k_loco[["1"]], addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],cores = 1)


pdf(file="~/Desktop/figs/4D.pdf", width = 10, height = 7)
plot_coefCC(ier5[5200:6700,], cross_basic$pmap,scan1_output = subset(x, lodcolumn="Ier5"), main = "Ier5 - Chr. 1", legend_ncol=1,top_panel_prop = 0.6)
plot_coefCC(qsox1[5200:6700,], cross_basic$pmap,scan1_output = subset(x, lodcolumn="Qsox1"), main = "Qsox1 - Chr. 1", legend_ncol=1,top_panel_prop = 0.6)
plot_coefCC(TMD_blup[5200:6700,], cross_basic$pmap,scan1_output = subset(DO_qtl_scan_normal, lodcolumn="uCT_Ct.TMD"), main = "TMD - Chr. 1", legend_ncol=1,top_panel_prop = 0.6)
plot_coefCC(ML_blup[5200:6700,], cross_basic$pmap,scan1_output = subset(DO_qtl_scan_normal, lodcolumn="ML"), main = "ML - Chr. 1", legend_ncol=1,top_panel_prop = 0.6)
dev.off()

##
##4E
#get microarray data from bioGPS (03.24.2020):
#GSE10246
biogps = read_csv("./data/GSE10246_bioGPS_032420.csv")
biogps = biogps[-which(as.numeric(biogps$`1420831_at`)<500),]
biogps$tis = rownames(biogps)
biogps = biogps[,-2]

#get mean of samples with more than one reading
u = unique(unlist(strsplit(biogps$tis, split = "[.]")))
u = u[-which(u %in% c("1","2"))]


biogps_means = as.data.frame(matrix(ncol=2, nrow=length(u)))
colnames(biogps_means) = c("val","tis")

for(i in 1:length(u)){
  sub = biogps[grep(x=biogps$tis, pattern = u[i]),]
  m = mean(as.numeric(sub[,1]))
  biogps_means[i,] = c(m, u[i])
}

biogps_means$col = 1
biogps_means[grep(pattern = "osteoblast", x = biogps_means$tis,ignore.case = T), "col"] = 2

p = biogps_means[order(as.numeric(biogps_means$val)),"tis"]

pdf(file="~/Desktop/figs/4D.pdf", width = 10, height = 7)
ggplot(biogps_means, aes(x=tis, y=as.numeric(val), fill=col)) + 
  geom_bar(position = "dodge",stat="identity") + coord_flip()  + xlab("Tissue") + ylab("Expression") + scale_color_brewer(palette = "Dark2") + theme(legend.position = "none") + scale_x_discrete(limits = p)
dev.off()
#
#
#


#6A
#calvarial ob data
#from B6_OBs_RNA_Seq.R





##6B
#seurat_analysis.R
load("./results/Rdata/seurat_ob.Rdata") #loads as "ob"

cairo_pdf(file="~/Desktop/figs/5Bqsox.pdf", width = 10, height = 7)
FeaturePlot(ob, features = "Qsox1", sort.cell = T, pt.size = 1)
dev.off()


cairo_pdf(file="~/Desktop/figs/5Bvim.pdf", width = 10, height = 7)
FeaturePlot(ob, features = "Vim", sort.cell = T, pt.size = 1)
dev.off()

cairo_pdf(file="~/Desktop/figs/5Blgal.pdf", width = 10, height = 7)
FeaturePlot(ob, features = "Lgals1", sort.cell = T, pt.size = 1)
dev.off()

cairo_pdf(file="~/Desktop/figs/5Bprrx.pdf", width = 10, height = 7)
FeaturePlot(ob, features = "Prrx2", sort.cell = T, pt.size = 1)
dev.off()


#6C


#Done externally











# #done externally
# library(BSgenome.Mmusculus.UCSC.mm10)
# from = 155778158
# to = 155812889
# txTr = GeneRegionTrack(TxDb.Mmusculus.UCSC.mm10.knownGene, chromosome = "chr1",start = from, end = to)
# plotTracks(txTr)
# 
# biomTrack <- BiomartGeneRegionTrack(genome = "mm10",name = "ENSEMBL", symbol = "Qsox1")
# plotTracks(biomTrack)
# 
# 
# gr=GRanges(seqnames = c("chr1"), ranges = IRanges(start=from, end=to), strand = strand(c("-")))
# subsetByOverlaps(transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene),gr)
# plot_theme = theme(legend.position = "none")
# geneViz(TxDb.Mmusculus.UCSC.mm10.knownGene, gr, BSgenome.Mmusculus.UCSC.mm10,isoformSel = "uc007dbp.1", labelTranscript = F, plotLayer = plot_theme)


##6D
#Mostly Charles Farber's code


dat2<-read.csv('./data/pheno_data/Qsox1_data/QsoxAssay_June09.csv',header=T)

dat2$pmol.H202.min.ul<-dat2$pmol.H2O2.synthesized...min/5

dat2$Genotype<-factor(dat2$Genotype, c('wt','Het','Mut'))
dat2<-filter(dat2,Genotype!='NA')
dat2.melt<-melt(dat2[,c(1,8,18)])
colnames(dat2.melt)<-c('Mutation','Genotype','variable','Activity')

cairo_pdf(file="~/Desktop/figs/6D.pdf", width = 10, height = 7)
ggplot(data = dat2.melt, aes(x=Genotype, y=Activity,Group=Mutation,fill=Mutation)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_point(  position="jitter", size=1) +
  ylab('QSOX1 Activity (pmol H2O2/min/ul')
dev.off()

####6E
# generate LSMEANS by sex adjusting for weight, length and mutation type
# can add mutation type to the filter statement to look at effects
# of individual mutation

dat1<-read.csv('./data/pheno_data/Qsox1_data/qsox1_caliper.csv',header=T,na.strings='NA')

dat1$Genotype<-factor(dat1$Genotype, levels=c('wt','Het','Mut'))
dat1<-dplyr::filter(dat1,dat1$Mom_genotype!='B6' & dat1$Dad_genotype!='B6')
dat1$ML_all<-(dat1$right_Femur_ML+dat1$left_Femur_ML)/2


dat1.m<-dplyr::filter(dat1,dat1$Sex=='M')


lf.lm.12<-lm(ML_all~Genotype+Weight+Qsox_Mutation+Length,data=dat1.m)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab='Femur Length (mm)')
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')


lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6Emales.pdf", width = 10, height = 7)
ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Het"),c("wt","Mut"),c("Het","Mut")), annotations=c("2.4e-3","2.87e-7","0.044"),aes(y=lsmeans.lsmean),y_position = c(1.95,1.98,1.97)) + ylim(1.83,1.98)+ylab("ML (mm)") + xlab("Genotype")
dev.off()

#females
dat1.f<-dplyr::filter(dat1,dat1$Sex=='F')


lf.lm.12<-lm(ML_all~Genotype+Weight+Qsox_Mutation+Length,data=dat1.f)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab='Femur Length (mm)')
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')


lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6E_females.pdf", width = 10, height = 7)
ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Het"),c("wt","Mut"),c("Het","Mut")), annotations=c("0.8","0.004","0.021"),aes(y=lsmeans.lsmean),y_position = c(1.70,1.74,1.73)) + ylim(1.64,1.74)+ylab("ML(mm)") + xlab("Genotype")
dev.off()

######
#ADD SAMPLE COUNTS PER GENOTYPE
######
#6F
#

dat1<-read.csv('./data/pheno_data/Qsox1_data/qsox_uCT.csv',header=T,na.strings='NA')

dat1$Genotype<-factor(dat1$Genotype, levels=c('wt','Het','Mut'))

lf.lm.12<-lm(pMOI.mm4.~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab=expression(paste("pMOI (",mm^4,")")))
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6Fpmoi.pdf", width = 10, height = 7)
ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.02"),aes(y=lsmeans.lsmean),y_position = c(0.52)) + ylim(0.36,0.52)+xlab("Genotype") + ylab("pMOI")
dev.off()
#
####
####

lf.lm.12<-lm(Imax.mm4.~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab=expression(paste("Imax (",mm^4,")")))
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6Fimax.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.009"),aes(y=lsmeans.lsmean),y_position = c(0.36)) + ylim(0.25,0.36)+xlab("Genotype") +ylab("Imax")
dev.off()
#
###
####

lf.lm.12<-lm(Ct.Ar.Tt.Ar....~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab=expression(paste("Ct.Ar.Tt.Ar.... (",mm^4,")")))
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6Fctarttar.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.031"),aes(y=lsmeans.lsmean),y_position = c(50)) + ylim(43,50)+xlab("Genotype") + ylab("Ct.Ar/Tt.Ar")
dev.off()
###
###
###

lf.lm.12<-lm(Tt.Ar.mm2.~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab=expression(paste("ttar (",mm^4,")")))
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6Fttar.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.083"),aes(y=lsmeans.lsmean),y_position = c(2.075)) + ylim(1.8,2.1)+xlab("Genotype") + ylab("Tt.Ar")
dev.off()
#
#
#
#

lf.lm.12<-lm(Ma.Ar.mm2.~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab="Ma.Ar")
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
pdf(file="~/Desktop/figs/6Fmaar.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.93"),aes(y=lsmeans.lsmean),y_position = c(1.09)) + ylim(0.98,1.09)+xlab("Genotype") + ylab("Ma.Ar")
dev.off()



#6G
lf.lm.12<-lm(Ct.TMD..mgHA.cm3.~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab="Ma.Ar")
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6G.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.4"),aes(y=lsmeans.lsmean),y_position = c(1182)) + ylim(1162,1182)+xlab("Genotype")+ylab("TMD")
dev.off()

#
#
#
###6H
lf.lm.12<-lm(Ct.Porosity....~Genotype+Weight..g.+Qsox.Mutation,data=dat1)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab="Ma.Ar")
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6H.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.242"),aes(y=lsmeans.lsmean),y_position = c(1.6)) + ylim(0.85,1.6)+xlab("Genotype") +ylab("Ct.Por")
dev.off()
#######
####6I
qsox_bend_AP = read.csv("./data/pheno_data/Qsox1_data/qsox1_bending_AP.csv")
qsox_bend_ML = read.csv("./data/pheno_data/Qsox1_data/qsox1_bending_ML.csv")

qsox_bend_ML$Genotype = qsox_bend_AP[match(qsox_bend_ML$Specimen, qsox_bend_AP$Specimen),"Genotype"]
qsox_bend_ML$Weight = qsox_bend_AP[match(qsox_bend_ML$Specimen, qsox_bend_AP$Specimen),"Weight"]
qsox_bend_ML$Qsox.Mutation = qsox_bend_AP[match(qsox_bend_ML$Specimen, qsox_bend_AP$Specimen),"Qsox.Mutation"]




###AP
lf.lm.12<-lm(Max.Load..N.~Genotype+Weight+Qsox.Mutation,data=qsox_bend_AP)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))

lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
cairo_pdf(file="~/Desktop/figs/6I_AP.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.198"),aes(y=lsmeans.lsmean),y_position = c(33)) + ylim(26,33)+xlab("Genotype") + ylab("Max Load (N)")
dev.off()

###

#6J
lf.lm.12<-lm(Max.Load..N.~Genotype+Weight+Qsox.Mutation,data=qsox_bend_ML)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
lff = lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')

lff = as.data.frame(lff)
pdf(file="~/Desktop/figs/6J_ML.pdf", width = 10, height = 7)

ggplot(lff,aes(x=lsmeans.Genotype)) + geom_errorbar(data = lff, aes(ymin=lsmeans.lower.CL, ymax=lsmeans.upper.CL),color="blue",size=3.5,width=0.15) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18))+geom_point(aes(y=lsmeans.lsmean),size=5)+geom_signif(comparisons = list(c("wt", "Mut")), annotations=c("0.001"),aes(y=lsmeans.lsmean),y_position = c(63)) + ylim(40,63)+xlab("Genotype") + ylab("Max Load (N)")

dev.off()


#6K is external







##SUPP FIG 1 in plot_overlaps_supplement.R

##SUPP fig 2
# #POMP data
# #from pomp_mapping_ML.R
load("./results/Rdata/scan1_pomp.Rdata")
load("./results/Rdata/coef_ML_pomp_blup.Rdata")
load("./data/POMP/MM_snps1_pmap.Rdata")

cairo_pdf(file="~/Desktop/figs/supp2_legend.pdf", width = 10, height = 7)
plot_coefCC(coef_ML_pomp_blup, MM_snps1_pmap["1"], scan1_output=subset(scan1_pomp, lodcolumn=2),legend = "topleft")
dev.off()
 

