# get data for one marker
read_intensities <-
    function(file="../RawData/intensities.fst", marker=NULL, snps=NULL)
{
    if(is.null(snps)) {
        snps <- fst::read.fst(file, column="snp")[,1]
    }

    if(is.null(marker)) { # if marker is NULL, return the snp names
        return(snps)
    } else {
        marker <- as.character(marker)
    }

    n <- match(marker, snps)
    if(is.na(n)) stop("marker ", marker, " not found")

    reorg_intensities( fst::read.fst(file, from=n, to=n+1) )
}


# reorganize SNP intensity data (for one marker)
reorg_intensities <-
    function(intensities, marker=NULL)
{
    if(!is.null(marker)) {
        marker <- as.character(marker)
        wh <- which(intensities[,1] == marker)
        if(sum(wh==0)) stop("Marker ", marker, " not found")
        intensities <- intensities[wh,]
    }

    result <- data.frame(X=as.numeric(intensities[1,-(1:2)]),
                         Y=as.numeric(intensities[2,-(1:2)]))
    rownames(result) <- colnames(intensities)[-(1:2)]

    result
}

# combine intensities and genotypes
grab_gni <-
    function(marker=NULL, cross=svenson, intensities_file="../RawData/intensities.fst",
             intensities_data=NULL, drop_bad_samples=TRUE)
{
    if(is.null(intensities_data)) {
        # get intensities
        int <- read_intensities(intensities_file, marker=marker)
    } else {
        int <- reorg_intensities(intensities_data, marker)
    }
    marker <- as.character(marker)


    # get genotypes
    chr <- qtl2::find_markerpos(cross, marker)$chr
    g <- cross$geno[[chr]][,marker]
    gg <- setNames(rep(0, nrow(int)), rownames(int))
    gg[names(g)] <- g

    if(drop_bad_samples) {
        int <- int[names(g),,drop=FALSE]
        gg <- gg[names(g)]
    }

    cbind(int, g=gg)
}

# load intensities and plot them, colored by genotype calls
#
# drop_bad_samples: if TRUE, don't plot points for samples that are not in the cross object
#                   (e.g., were omitted previously as being bad DNAs)
#
plot_intensities <-
    function(marker, cross=svenson, intensities_file="../RawData/intensities.fst",
             intensities_data=NULL, drop_bad_samples=TRUE, geno=NULL, ...)
{
    if(is.character(marker) || is.factor(marker)) {
        marker <- as.character(marker)
        gni <- grab_gni(marker, cross=cross, intensities_file=intensities_file,
                        intensities_data=intensities_data, drop_bad_samples=drop_bad_samples)
    } else {
        gni <- marker
    }

    if(!is.null(geno)) {
        if(is.logical(geno)) geno <- geno + 1 # FALSE/TRUE -> 1/2
        gni[names(geno),"g"] <- geno
    }

    internal_plot <-
        function(pch=21, bg=broman::brocolors("f2"),
                 xlab="allele 1", ylab="allele 2",
                 xlim=c(0, max(gni$X, na.rm=TRUE)),
                 ylim=c(0, max(gni$Y, na.rm=TRUE)),
                 ...)
        {
            grayplot(gni$X, gni$Y, pch=pch, bg=bg[match(gni$g, c(1:3,0))],
                     xlab=xlab, ylab=ylab,
                     xlim=xlim, ylim=ylim, ...)
        }

    internal_plot(...)

    invisible(gni)
}

# recode SNPs so that 1 = major allele in founders
recode_snps <-
    function(cross)
{
    for(i in seq_along(cross$geno)) {
        g <- cross$geno[[i]]
        g[g==0] <- NA
        fg <- cross$founder_geno[[i]]
        fg[fg==0] <- NA
        fgf <- colMeans(fg==3, na.rm=TRUE)
        if(any(!is.na(fgf) & fgf > 0.5)) {
            g[,fgf>0.5] <- 4 - g[,fgf>0.5]
            fg[,fgf>0.5] <- 4 - fg[,fgf>0.5]
        }
        g[is.na(g)] <- 0
        fg[is.na(fg)] <- 0
        cross$geno[[i]] <- g
        cross$founder_geno[[i]] <- fg
    }
    cross
}

# x,y -> sum,angle
transform_intensities <-
    function(intensities)
{
    rs <- sqrt(rowSums(intensities^2))
    cbind(sum=rs, angle=acos(intensities[,1]/rs))
}

test_founder_geno <-
    function(marker, cross=svenson, geno=m)
{
    chr <- find_markerpos(cross, marker)$chr
    fg <- cross$founder_geno[[chr]][,marker]

    cross <- cross[,chr]
    geno <- geno[,chr]
#    cross <- pull_markers(cross, marker)
#    geno[[chr]] <- geno[[chr]][,marker,drop=FALSE]

    nf <- length(fg)
    results <- vector("list", nf+1)

    results[[1]] <- qtl2::predict_snpgeno(cross, geno)[[chr]][,marker]

    for(i in 1:nf) {
        this_fg <- fg
        this_fg[i] <- 4-this_fg[i]
        cross$founder_geno[[chr]][,marker] <- this_fg
        results[[i+1]] <- qtl2::predict_snpgeno(cross, geno)[[chr]][,marker]
    }

    results <- matrix(unlist(results), ncol=nf+1)
    dimnames(results) <- list(rownames(geno[[1]]), c("-", names(fg)))
    results
}
