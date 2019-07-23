library(ggplot2)
library(qtl2)
##plots for paper##
load("./results/Rdata/apr_basic_cleaned.Rdata")
load("./results/Rdata/cross_basic_cleaned.Rdata")



#freq of founder genotype 
af_ind <- calc_geno_freq(apr, "individual")#marker

par(mar=c(4.1, 4.1, 0.6, 0.6))

par(mfrow=c(4,2))
for(i in 1:8){ hist(af_ind[,i], main=NULL, breaks=30,xlab = "geno freq by ind.",cex.lab=1.75, cex.axis=1.75)
  abline(v=mean(af_ind[,i]),col="red",lwd=3)
}

## allele freq per geno
#calculate allele freq of each marker, get chromosome from find_markerpos, plot contribution of allele per chrom
chr = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","X")
af_mar <- calc_geno_freq(apr, "marker", omit_x=FALSE)
af_mar = as.data.frame(af_mar)
af_mar$chrom = NA
af_mar$markerName = rownames(af_mar)


af_mar$chrom = as.character(apply(af_mar,1,function(x) find_markerpos(cross_basic$pmap,(x[10]))[1][[1]]))

#pseudos = grep(".loc",af_mar$markerName)
#for(i in pseudos){
#  af_mar[i,"chrom"] = strsplit(strsplit(af_mar[i,"markerName"],"\\.")[[1]][1],"c")[[1]][2]
#}

mean_per_chrom = aggregate(af_mar[,c(1:8)],by = list(af_mar$chrom), FUN = mean)

mean_per_chrom = mean_per_chrom[match(chr,mean_per_chrom$Group.1),]
rownames(mean_per_chrom) =mean_per_chrom$Group.1
mean_per_chrom = mean_per_chrom[,-1]

barplot(as.matrix(t(mean_per_chrom)),col = CCcolors, ylab = "Allele Frequency", xlab="Chr", main = "Global Allele Freq. Per Chromosome")
legend("topleft", fill=CCcolors, legend=c("A","B","C","D","E","F","G","H"))


##XO per generation
totxo <- as.data.frame(rowSums(nxo))
totxo$ngen = cross_basic$covar$ngen
totxo$ID = rownames(totxo)
totxo = totxo[order(as.numeric(totxo$ID)),]


g = ggplot(totxo, aes(x = seq_along(ID), y=`rowSums(nxo)`)) + geom_point(aes(color=ngen))+ geom_point(shape = 1, colour="black")
g + scale_colour_brewer(palette = "Paired")
#####
