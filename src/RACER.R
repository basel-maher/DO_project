#library(devtools)
#install_github("oliviasabik/RACER") 
library(RACER)
library(tidyverse)
library(data.table)

bmd = read.delim("~/Desktop/gefos_FNBMD_pooled_allChrom.txt",stringsAsFactors = FALSE)
lsbmd = read.delim("~/Desktop/gefos_LSBMD_pooled_allChrom.txt",stringsAsFactors = FALSE)

#gtex = read.delim("~/Desktop/Artery_Tibial", stringsAsFactors = F, header = FALSE, sep = " ")
#gtex = read.delim("~/Desktop/Skin_Sun_Exposed_Lower_leg", stringsAsFactors = F, header = FALSE, sep = " ")#800kb region missing
gtex = read.delim("~/Desktop/Skin_Sun_Exposed_Lower_leg83", stringsAsFactors = F, header = FALSE, sep = "\t")#800kb region missing
gtex = fread("~/Desktop/Heart_Left_Ventricle", stringsAsFactors = F, header = FALSE, sep = "\t",sep2 = " ")#800kb region missing


#gtex = read.delim("~/Desktop/GWAS_project/morris_ukbb/GTEx_processed/Skin_Sun_Exposed_Lower_leg", stringsAsFactors = F, header = FALSE, sep = " ")#800kb region missing

gtex_brain = read.delim("~/Desktop/Brain_Frontal_Cortex_BA9", stringsAsFactors = F, header = FALSE, sep = " ")

#gtex = gtex[which(gtex$V9 == "RHPN2"),]
gtex = gtex[which(gtex$V14 == "ACP2"),]
gtex = gtex[grep("ENSG00000134575", gtex$V1),]#ACP2

df <- gsub(pattern=" ", replacement="\t", gtex$V9, fixed = TRUE)
df <- do.call(rbind.data.frame, strsplit(df, split = "\t"))

gtex <- cbind(gtex[, c(1:8)],
                df,
                df[, c(1:ncol(df))])

lookup = read.delim("~/Desktop/chr11_lookup",stringsAsFactors = FALSE, header=F)

gtex$rsid = lookup[match(gtex$V2, lookup$V3),"V7"]
gtex$pos = lookup[match(gtex$V2, lookup$V3),"V2"]
gtex$chr = 11



gtex_brain = gtex_brain[which(gtex_brain$V9 == "SOST"),]

bmd = bmd[which(bmd$MarkerName %in% gtex$rsid),]

lsbmd = lsbmd[which(lsbmd$MarkerName %in% gtex$rsid),]

gtex_brain = gtex_brain[which(gtex_brain$V3 %in% lsbmd$MarkerName),]
gtex_brain = gtex_brain[-which(duplicated(gtex_brain$V3)),]
lsbmd$chr=17

gtex = gtex[which(gtex$rsid %in% bmd$MarkerName),]
#gtex = gtex[-which(duplicated(gtex$rsid)),]
bmd$chr = 11

bmd$pos = gtex[match(bmd$MarkerName,gtex$rsid),"pos"]

#for(i in 1:nrow(bmd)){
#  bmd$pos[i] = strsplit(bmd$pos[i], split = "_")[[1]][2]
#}

bmd_racer = RACER::formatRACER(assoc_data = bmd, chr_col = 8, pos_col = 9, p_col = 4)
#
#for(i in 1:nrow(gtex)){
#  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
#}

gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 12, pos_col = 11, p_col = 7)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 1, pops = "EUR", lead_snp = "rs7932354"))#rs10410778 for rhpn2

lead_snp = "rs7932354"


for (i in  1:nrow(bmd_racer_ld)) {
  tryCatch({

    if(is.na(bmd_racer_ld$LD[i])){
      ld_command = paste0('https://ldlink.nci.nih.gov/LDlinkRest/ldpair?var1=',lead_snp,'&var2=',bmd_racer_ld$RS_ID[i],'&pop=EUR&token=c0f613f149ab')
      z = as.data.frame(data.table::fread(ld_command,sep = "\t"))
      idx = grep("R2",z[,1])
      bmd_racer_ld$LD[i] = as.numeric(unlist(strsplit(z[idx,],split = " "))[2])
    }


  }, error=function(e){bmd_racer_ld$LD[i]=NA})
}

bmd_racer_ld$LD_BIN = cut(bmd_racer_ld$LD, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.2-0.0", "0.4-0.2", "0.6-0.4", "0.8-0.6", "1.0-0.8"))

# 
 gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 10, pops = "EUR", lead_snp = "rs7932354"))
# 
for (i in  1:nrow(gtex_racer_ld)) {
  tryCatch({

    if(is.na(gtex_racer_ld$LD[i])){
      ld_command = paste0('https://ldlink.nci.nih.gov/LDlinkRest/ldpair?var1=',lead_snp,'&var2=',bmd_racer_ld$RS_ID[i],'&pop=EUR&token=c0f613f149ab')
      z = as.data.frame(data.table::fread(ld_command,sep = "\t"))
      idx = grep("R2",z[,1])
      gtex_racer_ld$LD[i] = as.numeric(unlist(strsplit(z[idx,],split = " "))[2])
    }


  }, error=function(e){gtex_racer_ld$LD[i]=NA})
}

gtex_racer_ld$LD_BIN = cut(gtex_racer_ld$LD, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0.2-0.0", "0.4-0.2", "0.6-0.4", "0.8-0.6", "1.0-0.8"))

z = mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 11, plotby = "coord", name1 = "FNBMD",name2 = "ACP2 expression - SkinSunExp",start_plot = 46686260, end_plot = 46900000)
####
#400 434
gtex_brain = read.delim("~/Desktop/Brain_Frontal_Cortex_BA9", stringsAsFactors = F, header = FALSE, sep = " ")


gtex_brain = gtex_brain[which(gtex_brain$V9 == "SOST"),]

lsbmd = lsbmd[which(lsbmd$MarkerName %in% gtex_brain$V3),]

gtex_brain = gtex_brain[which(gtex_brain$V3 %in% lsbmd$MarkerName),]
gtex_brain = gtex_brain[-which(duplicated(gtex_brain$V3)),]
lsbmd$chr=17


lsbmd$pos = gtex_brain[match(lsbmd$MarkerName,gtex_brain$V3),"V2"]

for(i in 1:nrow(lsbmd)){
  lsbmd$pos[i] = strsplit(lsbmd$pos[i], split = "_")[[1]][2]
}

bmd_racer = RACER::formatRACER(assoc_data = lsbmd, chr_col = 8, pos_col = 9, p_col = 4)
#
for(i in 1:nrow(gtex_brain)){
  gtex_brain$pos[i] = strsplit(gtex_brain$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex_brain, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 1, pops = "EUR", lead_snp = "rs4792909"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs4792909"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 17, plotby = "gene", gene_plot = "SOST",name1 = "LSBMD", name2 = "SOST eQTL - Brain_Frontal_Cortex")






####

gtex_tibial = read.delim("~/Desktop/Artery_Tibial", stringsAsFactors = F, header = FALSE, sep = " ")
gtex_tibial = gtex_tibial[which(gtex_tibial$V9 == "TMEM263"),]

ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

ebmd = ebmd[which(ebmd$RSID %in% gtex_tibial$V3),]
gtex_tibial = gtex_tibial[which(gtex_tibial$V3 %in% ebmd$RSID),]

#gtex_tibial = gtex_tibial[-which(duplicated(gtex_tibial$V3)),]



bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 13)
#
for(i in 1:nrow(gtex_tibial)){
  gtex_tibial$pos[i] = strsplit(gtex_tibial$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex_tibial, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs4964511"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs4964511"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 12, plotby = "coord", gene_plot = "C12orf23",name1 = "eBMD", name2 = "TMEM263 eQTL - Artery_Tibial",start_plot = 107250000, end_plot = 107410000)
####















gtex_skin = read.delim("~/Desktop/Skin_Not_Sun_Exposed_Suprapubic", stringsAsFactors = F, header = FALSE, sep = " ")
gtex_skin = gtex_skin[which(gtex_skin$V9 == "JOSD1"),]

ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

ebmd = ebmd[which(ebmd$RSID %in% gtex_skin$V3),]
gtex_skin = gtex_skin[which(gtex_skin$V3 %in% ebmd$RSID),]

#gtex_skin = gtex_skin[-which(duplicated(gtex_skin$V3)),]



bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 13)
#
for(i in 1:nrow(gtex_skin)){
  gtex_skin$pos[i] = strsplit(gtex_skin$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex_skin, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs7290979"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs7290979"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 22, plotby = "coord", gene_plot = "JOSD1",name1 = "eBMD", name2 = "TMEM263 eQTL - Artery_Tibial",start_plot = 39000000, end_plot = 39200000,build = "hg19")
####
##
###
#
#
#
#
#
#
#
library(RACER)
#bmd = read.delim("~/Desktop/gefos_FNBMD_pooled_allChrom.txt",stringsAsFactors = FALSE)
#lsbmd = read.delim("~/Desktop/gefos_LSBMD_pooled_allChrom.txt",stringsAsFactors = FALSE)
ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("~/Desktop/Adipose_Subcutaneous", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V9 == "SNN"),]

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]


gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]
gtex = gtex[-which(duplicated(gtex$V3)),]



ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"V2"]


bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 13)
bmd_racer = bmd_racer[,-11]
#
for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs7102"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs7102"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 16, plotby = "gene",gene_plot = "SNN",build = "hg19", name1 = "eBMD", name2 = "Snn - Adipose Subcutaneous")
####





ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("~/Desktop/Artery_Tibial", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V9 == "ITIH4"),]

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]


gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]
#gtex = gtex[-which(duplicated(gtex$V3)),]



ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"V2"]


bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 13)
bmd_racer = bmd_racer[,-11]
#
for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs13072536"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs13072536"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 3, plotby = "gene",gene_plot = "ITIH4",build = "hg19", name1 = "eBMD", name2 = "Itih4 - Tibial Artery")
####






ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("~/Desktop/Artery_Tibial", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V9 == "RHPN2"),]

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]


gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]
#gtex = gtex[-which(duplicated(gtex$V3)),]



ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"V2"]

ebmd = ebmd[,-11]

ebmd = unique(ebmd)
gtex = unique(gtex)
bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 12)
bmd_racer = bmd_racer[,-14]
#
for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs11881367"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs11881367"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 19, plotby = "gene",gene_plot = "RHPN", build = "hg19", name1 = "eBMD", name2 = "Acp2 - Skin_Sun_Exposed")
####











ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("~/Desktop/Heart_Left_Ventricle", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V9 == "FAT1"),]

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]


gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]
#gtex = gtex[-which(duplicated(gtex$V3)),]



ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"V2"]

ebmd = ebmd[,-11]

ebmd = unique(ebmd)
gtex = unique(gtex)
bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 12)
bmd_racer = bmd_racer[,-14]
#
for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs327102"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs327102"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 4,plotby = "coord", build = "hg19", name1 = "eBMD", name2 = "FAT1 - Heart_Left_Ventricle", start_plot = 187500000, end_plot = 187800000)
####









ebmd = read.delim("~/Desktop/GWAS_project/morris_ukbb/Morrisetal2018.NatGen.SumStats/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt",stringsAsFactors = FALSE)

gtex = read.delim("~/Desktop/Brain_Nucleus_accumbens_basal_ganglia", stringsAsFactors = F, header = FALSE, sep = " ")

gtex = gtex[which(gtex$V9 == "CSNK1G3"),]

ebmd = ebmd[which(ebmd$RSID %in% gtex$V3),]


gtex = gtex[which(gtex$V3 %in% ebmd$RSID),]
#gtex = gtex[-which(duplicated(gtex$V3)),]



ebmd$pos = gtex[match(ebmd$RSID,gtex$V3),"V2"]

ebmd = ebmd[,-11]

ebmd = unique(ebmd)
gtex = unique(gtex)
bmd_racer = RACER::formatRACER(assoc_data = ebmd, chr_col = 3, pos_col = 4, p_col = 12)
bmd_racer = bmd_racer[,-14]
#
for(i in 1:nrow(gtex)){
  gtex$pos[i] = strsplit(gtex$V2[i], split = "_")[[1]][2]
}

gtex_racer = RACER::formatRACER(assoc_data = gtex, chr_col = 11, pos_col = 13, p_col = 5)

bmd_racer_ld = (RACER::ldRACER(assoc_data = bmd_racer, rs_col = 2, pops = "EUR", lead_snp = "rs7703751"))
gtex_racer_ld = (RACER::ldRACER(assoc_data = gtex_racer, rs_col = 3, pops = "EUR", lead_snp = "rs7703751"))


mirrorPlotRACER(assoc_data1 = bmd_racer_ld, assoc_data2 = gtex_racer_ld, chr = 5,plotby = "coord", build = "hg19", name1 = "estimated BMD", name2 = "CSNK1G3 - Brain_Nuc_Acc_Basal_Ganglia",start_plot = 122650000, end_plot = 123050000)
####
