#FROM https://github.com/kbroman/Paper_MPPdiag/edit/master/R/geneseek2fst.R

# grab intensities from GeneSeek FinalReport.txt files
# convert them to a single big data frame, and save for use with fst

library(fst)

# simple version of data.table::fread()
myfread <- function(filename) data.table::fread(filename, data.table=FALSE, skip=9)

ifiles <- c("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20171110_FinalReport.txt",
            "./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20180408_FinalReport.txt",
            "./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20181213_FinalReport.txt",
            "./data/GIGAMUGA/FinalReport_files/20190228_FinalReport_FIXED.txt")

dat <- NULL
for(file in ifiles) {
  if(is.null(dat)) dat <- myfread(file)[,c("SNP Name", "Sample ID", "X", "Y")]
  else dat <- rbind(dat, myfread(file)[,c("SNP Name", "Sample ID", "X", "Y")])
}

# create matrices that are snps x samples
snps <- unique(dat[,1])
samples <- unique(dat[,2])
X <- Y <- matrix(ncol=length(samples), nrow=length(snps))
dimnames(X) <- dimnames(Y) <- list(snps, samples)
for(i in seq(along=samples)) {
  cat(i, "of", length(samples), "\n")
  tmp <- dat[dat[,2]==samples[i],]
  X[,as.character(samples[i])] <- tmp[,3]
  Y[,as.character(samples[i])] <- tmp[,4]
}

# bring together in one matrix
result <- cbind(snp=rep(snps, 2),
                channel=rep(c("X", "Y"), each=length(snps)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps), seq_along(snps)+length(snps)))),]
rownames(result) <- 1:nrow(result)

# write to fst file
write.fst(result, "./results/GIGAMUGA/intensities.fst", compress=100)
