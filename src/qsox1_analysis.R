library(lsmeans)



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

