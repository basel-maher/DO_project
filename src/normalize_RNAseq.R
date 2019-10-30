#This script normalizes RNA-seq counts#
library(DESeq2)


#read in the RNA-seq counts
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#variance-stabilize the counts
dds.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))

#add ".1" to some colnames so they match the GigaMUGA data
# samples 4 and 371 need ".1" added
colnames(dds.vst)[which(colnames(dds.vst) == 4)] = "4.1"
colnames(dds.vst)[which(colnames(dds.vst) == 371)] = "371.1"

#transpose
dds.vst = t(dds.vst)
dds.vst_n = as.data.frame(dds.vst)

#perform quantile-based inverse normal transform. aka match each gene datapoint to a quantile, then match to quantile in normal distribution
#from https://www.nature.com/articles/nature11401#s1 (FTO genotype BMI Visscher et al)

for(col in 1:ncol(dds.vst)){
  dds.vst[,col] = qnorm((rank(dds.vst[,col],na.last="keep")-0.5)/sum(!is.na(dds.vst[,col])))
}

dds.vst = as.data.frame(dds.vst)

dds.vst$Mouse.ID = rownames(dds.vst)
dds.vst = dds.vst[,c(ncol(dds.vst), 1:(ncol(dds.vst)-1))]
#write
write.table(dds.vst, "./results/flat/counts_vst_qnorm.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
save(object = dds.vst, file = "./results/Rdata/counts_vst_qnorm.Rdata")

dds.vst = dds.vst[,-1]
write.table(dds.vst, "./results/flat/counts_vst_qnorm_nohead.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")


