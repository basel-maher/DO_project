
# read in data and order genotype classifications (this will make sure
# genotypes are listed in order)
dat1<-read.csv('./data/pheno_data/Qsox1_data/Qsox_31July2018.csv',header=T,na.strings='NA')
dat1$Genotype<-factor(dat1$Genotype, levels=c('wt','Het','Mut'))

# count number of mice by genotype, sex and mutation type
table(dat1$Qsox_Mutation,dat1$Genotype,dat1$Sex)

# create derived variables
library(dplyr)
dat1$left_MLAP_ratio<-dat1$left_Femur_ML/dat1$left_Femur_AP
dat1$right_MLAP_ratio<-dat1$right_Femur_ML/dat1$right_Femur_AP
dat1$left_MLxAP<-dat1$left_Femur_ML+dat1$left_Femur_AP
dat1$right_MLxAP<-dat1$right_Femur_ML+dat1$right_Femur_AP
dat1$ell<-sqrt((dat1$right_Femur_ML^2-dat1$right_Femur_AP^2)/dat1$right_Femur_ML^2)
dat1$ML_all<-(dat1$right_Femur_ML+dat1$left_Femur_ML)/2
dat1$FL_all<-(dat1$right_Femur_FL+dat1$left_Femur_FL)/2
dat1$AP_all<-(dat1$right_Femur_AP+dat1$left_Femur_AP)/2
dat1$Tt.Ar<-((dat1$FL_all/2)*(dat1$AP_all/2)*3.1415)


# QC all phenotypes with histograms
hist(dat1$ML_all)

# eliminate samples from F1 mice
dat1<-filter(dat1,dat1$Mom_genotype!='B6' & dat1$Dad_genotype!='B6')

# run ANOVAs, adjusting for weight and mutation type
anova(lm((ML_all)~Genotype+Weight+Sex,data=dat1))

# generate means
tapply(dat1$ML_all,dat1$Genotype,mean,na.rm=T)
boxplot(dat1$ML_all~dat1$Genotype,ylab='Femur Length (mm)',
        xlab='Qsox1 Genotype',main='P=0.0516')

# generate LSMEANS by sex adjusting for weight and mutation type
# can add mutation type to the filter statement to look at effects
# of individual mutation
library(lsmeans)
library(ggplot2)
rm(dat1.s)
table(dat1$Qsox_Mutation)
dat1.s<-filter(dat1,dat1$Sex=='M' & (dat1$Qsox_Mutation=='del_1' | dat1$Qsox_Mutation=='del_7+6' |
                                       dat1$Qsox_Mutation=='del_171' | dat1$Qsox_Mutation=='del_171_SJL' | 
                                       dat1$Qsox_Mutation=='del_1300' | dat1$Qsox_Mutation=='del_1300_SJL' |
                                       dat1$Qsox_Mutation=='del_780' ))

tapply(dat1.s$Length,dat1.s$Genotype,mean,na.rm=T)

table(dat1.s$Genotype);sum(table(dat1.s$Genotype))
table(dat1.s$Genotype,dat1.s$Qsox_Mutation)

ggplot(dat1.s,aes(x=Genotype, y=Length)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.2))
hist(dat1.s$Length)
table(dat1.s$Genotype)
colnames(dat1.s)
lf.lm.12<-lm(ML_all~Genotype+Weight+Qsox_Mutation,data=dat1.s)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab='Femur Length (mm)')
lsmeans(lf.lm.12,pairwise~Genotype, adj='none')

# power calculations
rm(p.dist.20)
p.dist.20<-c()
for(i in 1:1000){
  dat1.s1<-sample(filter(dat1.s,Genotype=='wt')$ML_all,25)
  dat1.s2<-sample(filter(dat1.s,Genotype=='Mut')$ML_all,25)
  dat1.s3<-rbind(dat1.s1,dat1.s2)
  p.dist.20<-c(p.dist.20,t.test(dat1.s1,dat1.s2)$p.value)
  print(i)
}
  hist(-log10(p.dist.20))
(length(p.dist.20[p.dist.20>=0.05])/1000)*100


# test of sampling
dat1.s<-filter(dat1,dat1$Sex=='F' & (dat1$Qsox_Mutation=='del_7+6' | dat1$Qsox_Mutation=='del_1bp' |
                                       dat1$Qsox_Mutation=='del_171' | dat1$Qsox_Mutation=='del_171_SJL' | dat1$Qsox_Mutation=='del_1300' |
                                       dat1$Qsox_Mutation=='del_780'))
dim(dat1.s)
ran.mean<-c()
for(i in 1:1000){
wt.s<-filter(dat1.s,dat1.s$Sex=='M' & dat1.s$Genotype=='wt')
het.s<-filter(dat1.s,dat1.s$Sex=='M' & dat1.s$Genotype=='Het')
mut.s<-filter(dat1.s,dat1.s$Sex=='M' & dat1.s$Genotype=='Mut')
wt.s<-wt.s[sample(c(1:nrow(wt.s)),30,replace=F),]
het.s<-het.s[sample(c(1:nrow(het.s)),30,replace=F),]
mut.s<-mut.s[sample(c(1:nrow(mut.s)),30,replace=F),]
all.s<-rbind(wt.s,het.s,mut.s)
z1<-tapply(all.s$ML_all,all.s$Genotype,mean,na.rm=T)
ran.mean[i]<-z1[3]-z1[1]
print(i)
}

##### plot Qsox1 assay data

dat2<-read.csv('./data/pheno_data/Qsox1_data/QsoxAssay_June09.csv',header=T)

head(dat2)
dat2$pmol.H202.min.ul<-dat2$pmol.H2O2.synthesized...min/5

dat2.f<-filter(dat2,Sex=='F')
tapply(dat2.f$pmol.H202.min.ul,dat2.f$Genotype,mean)

dat2.m<-filter(dat2,Sex=='M')
tapply(dat2.m$pmol.H202.min.ul,dat2.m$Genotype,mean)

library(reshape2)
colnames(dat2)
dat2$Genotype<-factor(dat2$Genotype, c('wt','Het','Mut'))
dat2<-filter(dat2,Genotype!='NA')
dat2.melt<-melt(dat2[,c(1,8,18)])
head(dat2.melt)
colnames(dat2.melt)<-c('Mutation','Genotype','variable','Activity')

ggplot(data = dat2.melt, aes(x=Genotype, y=Activity,Group=Mutation,fill=Mutation)) + 
  geom_boxplot(outlier.shape=NA) +
  geom_point(  position="jitter", size=1) +
  ylab('QSOX1 Activity (pmol H2O2/min/ul')


p <- ggplot(mtcars, aes(factor(cyl), mpg)) + 
  geom_boxplot(outlier.shape=10, outlier.size=8)  +
  geom_point(aes(factor(cyl), mpg, color=mpg),  position="jitter", size=1)
p



#######################
qsox1_uct = read.csv('./data/pheno_data/Qsox1_data/qsox_uCT.csv',header=T,na.strings='NA')

lf.lm.12<-lm(Ct.Porosity....~Genotype+Weight..g.+Qsox.Mutation,data=qsox1_uct)
anova(lf.lm.12)
lf<-lsmeans(lf.lm.12,"Genotype")
par(mfrow=c(2,2))
plot(lf,horizontal=F,ylab='Femur Length (mm)')
lsmeans(lf.lm.12,pairwise~Genotype, adj='tukey')



#ct.ar, ct.ar/tt.ar, ct.th, pMOI, Imax, 


