#figures/tables for paper

set.seed(8675309)
library(qtl2)
library(ggplot2)
library(gridExtra)
library(grid)

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

par(mar=c(4.1, 4.1, 0.6, 0.6))

par(mfrow=c(4,2))
for(i in 1:8){ hist(af_ind[,i], main=NULL, breaks=30,xlab = "geno freq by ind.",cex.lab=1.25, cex.axis=1.25)
  abline(v=mean(af_ind[,i]),col="red",lwd=3)
}



#1C)
#boxplot BV/TV
#make dataframe with bv/tv and sex
df_bvtv = as.data.frame(cross_basic$pheno[,"uCT_BV.TV"])
df_bvtv$pheno = "BV/TV"
df_bvtv$ind = rownames(df_bvtv)
df_bvtv$M = covar[,"sex"]
colnames(df_bvtv)[1] = "val"

df_str = as.data.frame(cross_basic$pheno[,"bending_max_load"])
df_str$pheno = "max_load"
df_str$ind = rownames(df_str)
df_str$M = covar[,"sex"]
colnames(df_str)[1] = "val"


df = rbind(df_bvtv, df_str)
df_2 = df

df$pheno = paste0(df$pheno,"_",df$M)
df = rbind(df,df_2)

p <- ggplot(df, aes(x=pheno, y=val)) + 
  geom_boxplot()

p
#1D) heritability

#all samples
hsq = est_herit(pheno = cross_basic$pheno, kinship = k, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")], max_batch=1)
hsq_df = as.data.frame(hsq)
#males
m = which(covar[,"sex"] == 1)
hsq_m = est_herit(pheno = cross_basic$pheno[m,], kinship = k, addcovar = covar[,c("age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")], max_batch=1)
hsq_m_df = as.data.frame(hsq_m)

#females
f = which(covar[,"sex"] == 0)
hsq_f = est_herit(pheno = cross_basic$pheno[f,], kinship = k, addcovar = covar[,c("age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")], max_batch=1)
hsq_f_df = as.data.frame(hsq_f)


hsq_table= merge(hsq_df, hsq_m_df, by=0)
rownames(hsq_table) = hsq_table$Row.names
hsq_table = hsq_table[,-1]

hsq_table= merge(hsq_table, hsq_f_df, by=0)

write.table(hsq_table, file="~/Desktop/hsq.csv", quote=F, sep=",")


grid.table(hsq_table)
hsq_table = hsq_table[-which(hsq_table$Row.names %in% c("sex", "sac_time", "age_at_sac_days","Mouse.ID")),]

hsq_table = hsq_table[order(hsq_table$hsq),]

hsq_table$Row.names <- factor(hsq_table$Row.names, levels = hsq_table$Row.names)

p<-ggplot(hsq_table, aes(x=Row.names, y=hsq)) + 
  geom_point() + coord_flip()
p

x = as.data.frame(hsq)
x$pheno = rownames(x)
x$sex = "B"

y = hsq_m_df
y$pheno = rownames(y)
y$sex = "M"
colnames(y)[1] = "hsq"

z = hsq_f_df
z$pheno = rownames(z)
z$sex = "F"
colnames(z)[1] = "hsq"

x = rbind(x,y)
x = rbind(x,z)
x = x[-which(x$pheno %in% c("sex", "sac_time", "age_at_sac_days","Mouse.ID")),]


p<-ggplot(x, aes(x=pheno, y=hsq)) + 
  geom_point(aes(color=sex, shape=sex)) + coord_flip()
#####
#####
#####


