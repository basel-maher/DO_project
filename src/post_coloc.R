# Post-coloc analyses

#homs = c(homologs_win_1_estrada, homologs_win_1_morris)
#homs=(unique(homs))

#544 overall BANs overlap gwas loci

##colocalizing eqtl
#morris_coloc = read.table("~/Desktop/coloc_BAN_over75.txt")
#estrada_coloc = read.table("~/Desktop/coloc_v7_all_results.txt")
#estrada_coloc = estrada_coloc[which(as.numeric(estrada_coloc$H4) >= 0.75),]

#coloc = c(morris_coloc$gene,estrada_coloc$gene)
#unique(coloc)

#MYPOP,PLEKHM1,ZNF609,GPR133,PTRF

#plekhm1 in superset. so is ptrf as cavin1. so is gpr133 as adgrd1

#q=32 #number of bone genes in coloc genes
#k=51 #number of coloc genes
#m=length(which(tolower(homs) %in% superset))#number bone genes in BAN homologs
#n=length(which(tolower(homs) %in% superset == FALSE))   #number of non-bone genes in BAN homologs

phyper(q,m,n,k,lower.tail = F)