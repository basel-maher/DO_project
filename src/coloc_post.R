# Post-coloc analyses

#How many BANs within 1 Mbp of lead snps?
hom_morris = read.table("./results/flat/homologs_within1mbp_morris_REV2.txt")
hom_estrada = read.table("./results/flat/homologs_within1_estrada_REV2.txt")

homologs = c(hom_estrada$V1, hom_morris$V1)
homologs = unique(homologs)

length(homologs)#688 genes

##How many lead snps?
estrada_lead_snps = read.table("./data/GEFOS/lead_snps_pos", header=F, stringsAsFactors = F)
morris_lead_snps = read.csv("./data/Morrisetal2018.NatGen.SumStats/Morris_eBMD_conditionally_ind_snps.csv", header=T, stringsAsFactors = F)

lead = c(estrada_lead_snps$V6, morris_lead_snps$RSID)
lead = unique(lead)
length(lead)#1161 unique lead snps


#Of the 688 homologs within 1 mbp of a lead snp, how many have colocalizing eQTL?
fn = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_FNBMD_over75_REV2.txt")
ls = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_LSBMD_over75_REV2.txt")
morris_coloc = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_morris_v7_all_results_over75_REV2.txt")

coloc_genes = c(fn$gene, ls$gene, morris_coloc$gene)
coloc_genes = unique(coloc_genes)
length(coloc_genes)#66


#How many are known regulators of bone biology?
superset = read.table("./results/flat/superset_in_networks.txt")
superset = superset$V1
length(which(tolower(coloc_genes) %in% tolower(superset)))
#34 in superset. however, GPR133 amd ZNF609 are also in superset as Adgrd1 and Zfp609, respectively. So 36 in superset


#based on overlap with bone list, are the colocalizing genes enriched in bone genes? Can do this relative to enrichment in genome or enrichment in BANs that were colocalized 

#This analysis considers the BAN genes as the global set
new_set = superset
new_set[which(tolower(new_set) == "cavin1")] = "ptrf"
new_set[which(tolower(new_set) == "adgrd1")] = "gpr133"
new_set[which(tolower(new_set) == "zfp609")] = "znf609"

length(which(tolower(coloc_genes) %in% tolower(new_set)))

allgenes = homologs

a = length(which(tolower(coloc_genes) %in% tolower(new_set))) # number of coloc genes that are also bone genes (36)
b = length(which(tolower(coloc_genes) %in% tolower(new_set) == F)) # number of coloc genes that are not bone genes (30, 66-a)
c = length(which((tolower(allgenes) %in% tolower(new_set)) & (tolower(allgenes) %in% tolower(coloc_genes) == FALSE)))#not coloc and bone genes (200)
d = length(which((tolower(allgenes) %in% tolower(new_set)==FALSE) & (tolower(allgenes) %in% tolower(coloc_genes) == FALSE)))#not coloc and not bone genes (422)

#make contingency table
mat = matrix(nrow=2,ncol=2)
mat[1,] = c(a,c)
mat[2,] = c(b,d)
fisher.test(mat) #two sided (enrichment or depletion)
fisher.test(mat,alternative = "g") #enrichment
(x=fisher.test(mat,alternative = "g")$p.value)
###same but with hypergeometric 


# #This analysis considers all genes in genome as the global set
# #23648 genes
# 
# new_set = superset
# new_set[which(tolower(new_set) == "cavin1")] = "ptrf"
# new_set[which(tolower(new_set) == "adgrd1")] = "gpr133"
# new_set[which(tolower(new_set) == "zfp609")] = "znf609"
# 
# 
# allgenes = homologs
# 
# a = length(which(tolower(coloc_genes) %in% tolower(new_set))) # number of coloc genes that are also bone genes (40)
# b = length(which(tolower(coloc_genes) %in% tolower(new_set) == F)) # number of coloc genes that are not bone genes (32, 72-a)
# c = length(which(tolower(new_set) %in% tolower(coloc_genes) == FALSE))#not coloc and bone genes. For the genome it would equal number of bone genes that dont colocalize (1251)
# d = 23645 - (a+b+c)#not coloc and not bone genes (22087)
# 
# #make contingency table
# mat = matrix(nrow=2,ncol=2)
# mat[1,] = c(a,c)
# mat[2,] = c(b,d)
# fisher.test(mat) #two sided
# fisher.test(mat,alternative = "g") #enrichment


