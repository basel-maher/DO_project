library(qtl2convert)
#ADAPTED FROM Preparing Diversity Outbred (DO) mouse data for R/qtl2 
#https://kbroman.org/qtl2/pages/prep_do_data.html

#This is adapted from geneseek2qtl2.R from https://kbroman.org/qtl2/assets/geneseek2qtl2.R

# convert GeneSeek FinalReport files to format for R/qtl2
#
# - creates one genotype CSV file for each chromosome
#
# - also creates 4 files containing the two channels of SNP intensities for markers on the X and Y chr
#   (these are useful for verifying the sex of the mice)

# file containing allele codes for GigaMUGA data
#   - from GM_processed_files.zip, https://doi.org/10.6084/m9.figshare.5404759

codefile <- "./data/GIGAMUGA/GM_processed_files/GM_allelecodes.csv"

# input files with GigaMUGA genotypes
#  - can be a single file or a vector of multiple files
#  - if samples appear in multiple files, the genotypes in later files
#    will be used in place of genotypes in earlier files
#  - files can be gzipped (".gz" extension)

#These are GeneSeek FinalReport files
ifiles <- c("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20171110_FinalReport.txt",
            "./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20180408_FinalReport.txt",
            "./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20181213_FinalReport.txt",
            "./data/GIGAMUGA/FinalReport_files/20190228_FinalReport_FIXED.txt")

# file "stem" for output files
# output files will be like "DO_qtl2_geno19.csv"
ostem <- "./results/GIGAMUGA/qtl2_batches_1-4/DO_qtl2"

##############################
# define a couple of functions
##############################
# simple version of data.table::fread()

#skip 9 to skip the header lines in FinalReports
myfread_fr <- function(filename) data.table::fread(filename, data.table=FALSE,skip = 9)
#skip 3 header lines in codefile
myfread_a <- function(filename) data.table::fread(filename, data.table=FALSE,skip = 3)


# cbind, replacing matching columns with second set and adding unique ones
cbind_smother <-
  function(mat1, mat2)
  {
    cn1 <- colnames(mat1)
    cn2 <- colnames(mat2)
    m <- (cn2 %in% cn1)
    if(any(m)) {
      mat1[,cn2[m]] <- mat2[,cn2[m],drop=FALSE]
      if(any(!m)) {
        mat1 <- cbind(mat1, mat2[,cn2[!m],drop=FALSE])
      }
    }
    else {
      mat1 <- cbind(mat1, mat2)
    }
    
    mat1
  }
##############################



# read genotype codes
#skip header
codes <- myfread_a(codefile)

full_geno <- NULL
cXint <- cYint <- NULL

for(ifile in ifiles) {
  cat(" -File:", ifile, "\n")
  rezip <- FALSE
  if(!file.exists(ifile)) {
    cat(" -Unzipping file\n")
    system(paste("gunzip", ifile))
    rezip <- TRUE
  }
  
  cat(" -Reading data\n")
  g <- myfread_fr(ifile)
  # subset to the markers in the codes object
  g <- g[g[,"SNP Name"] %in% codes[,"marker"],]
  
  # NOTE: may need to revise the IDs in the 2nd column
  samples <- unique(g[,"Sample ID"])
  
  # matrix to contain the genotypes
  geno <- matrix(nrow=nrow(codes), ncol=length(samples))
  dimnames(geno) <- list(codes[,"marker"], samples)
  
  # fill in matrix
  cat(" -Reorganizing data\n")
  for(i in seq(along=samples)) {
    if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
    wh <- (g[,"Sample ID"]==samples[i])
    geno[g[wh,"SNP Name"],i] <- paste0(g[wh,"Allele1 - Forward"], g[wh,"Allele2 - Forward"])
  }
  
  cat(" -Encode genotypes\n")
  geno <- qtl2convert::encode_geno(geno, as.matrix(codes[,c("A","B")]))
  
  if(is.null(full_geno)) {
    full_geno <- geno
  } else {
    # if any columns in both, use those from second set
    #full_geno <- cbind_smother(full_geno, geno)
    cn1 <- colnames(full_geno)
    cn2 <- colnames(geno)
    m <- (cn2 %in% cn1)
    if(any(m)) {
      full_geno[,cn2[m]] <- geno[,cn2[m],drop=FALSE]
      if(any(!m)) {
        full_geno <- cbind(full_geno, geno[,cn2[!m],drop=FALSE])
      }
    }
    else {
      full_geno <- cbind(full_geno, geno)
    }
  }
  
  # grab X and Y intensities
  cat(" -Grab X and Y intensities\n")
  gX <- g[g[,"SNP Name"] %in% codes[codes$chr=="X","marker"],]
  gY <- g[g[,"SNP Name"] %in% codes[codes$chr=="Y","marker"],]
  cX <- matrix(nrow=sum(codes$chr=="X"),
               ncol=length(samples))
  dimnames(cX) <- list(codes[codes$chr=="X","marker"], samples)
  cY <- matrix(nrow=sum(codes$chr=="Y"),
               ncol=length(samples))
  dimnames(cY) <- list(codes[codes$chr=="Y","marker"], samples)
  for(i in seq(along=samples)) {
    if(i==round(i,-1)) cat(" --Sample", i, "of", length(samples), "\n")
    wh <- (gX[,"Sample ID"]==samples[i])
    cX[gX[wh,"SNP Name"],i] <- (gX$X[wh] + gX$Y[wh])/2
    
    wh <- (gY[,"Sample ID"]==samples[i])
    cY[gY[wh,"SNP Name"],i] <- (gY$X[wh] + gY$Y[wh])/2
  }
  if(is.null(cXint)) {
    cXint <- cX
    cYint <- cY
  } else {
    # if any columns in both, use those from second set
    #cXint <- cbind_smother(cXint, cX)
    cn1 <- colnames(cXint)
    cn2 <- colnames(cX)
    m <- (cn2 %in% cn1)
    if(any(m)) {
      cXint[,cn2[m]] <- cX[,cn2[m],drop=FALSE]
      if(any(!m)) {
        cXint <- cbind(cXint, cX[,cn2[!m],drop=FALSE])
      }
    }
    else {
      cXint <- cbind(cXint, cX)
    }
    #cYint <- cbind_smother(cYint, cY)
    cn1 <- colnames(cYint)
    cn2 <- colnames(cY)
    m <- (cn2 %in% cn1)
    if(any(m)) {
      cYint[,cn2[m]] <- cY[,cn2[m],drop=FALSE]
      if(any(!m)) {
        cYint <- cbind(cYint, cY[,cn2[!m],drop=FALSE])
      }
    }
    else {
      cYint <- cbind(cYint, cY)
    }
  }
  
  if(rezip) {
    cat(" -Rezipping file\n")
    system(paste("gzip", ifile))
  }
}

# write X and Y intensities
cat(" -Writing X and Y intensities\n")
qtl2convert::write2csv(cbind(marker=rownames(cXint), cXint),
                       paste0(ostem, "_chrXint.csv"),
                       paste(ostem, "X chr intensities"),
                       overwrite=TRUE)
qtl2convert::write2csv(cbind(marker=rownames(cYint), cYint),
                       paste0(ostem, "_chrYint.csv"),
                       paste(ostem, "Y chr intensities"),
                       overwrite=TRUE)

# write data to chromosome-specific files
cat(" -Writing genotypes\n")
for(chr in c(1:19,"X","Y","M")) {
  mar <- codes[codes$chr==chr,"marker"]
  g <- full_geno[mar,]
  qtl2convert::write2csv(cbind(marker=rownames(g), g),
                         paste0(ostem, "_geno", chr, ".csv"),
                         paste0(ostem, " genotypes for chr ", chr),
                         overwrite=TRUE)
}

