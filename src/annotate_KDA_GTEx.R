options(stringsAsFactors = FALSE)


#function to convert mouse to human gene names
#use mgi homolog list
#http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt (Mar-14-2020)

homology = read.table("./data/mgi_homologs.txt",sep = "\t", header = T)

convertMousetoHuman = function(x){
  human=c()
  for(i in 1:length(x)){
    id = homology[which((homology$Common.Organism.Name == "mouse, laboratory") & tolower(homology$Symbol) == tolower(x[i])),"HomoloGene.ID"]
    hum = homology[which((homology$Common.Organism.Name == "human") & homology$HomoloGene.ID == id),"Symbol"]
    human=append(human,hum)
  }
  human=unique(human)
  return(human)
 
}






#Get relevant genes from kda
kda = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv")
kda_m=read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv")
kda_f = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv")

nominal = kda[which(kda$hyper <= 0.05),"gene"]
nominal = append(nominal,kda_m[which(kda_m$hyper <= 0.05),"gene"])
nominal = append(nominal,kda_f[which(kda_f$hyper <= 0.05),"gene"])

nominal=unique(nominal) #1050 unique nominal key drivers

#convert mouse genes to human equivalent

nominal_hum = convertMousetoHuman(nominal) #904


##get GTEx v7 output files

morris_75 = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_morris_v7_all_results_over75.txt")
fnbmd_75 = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_FNBMD_over75.txt")
lsbmd_75 = read.table("~/Documents/projects/DO_project/results/flat/coloc/coloc_v7_LSBMD_over75.txt")


all_coloc_greaterorover75 = rbind(morris_75,fnbmd_75,lsbmd_75)

write.csv(all_coloc_greaterorover75,file = "~/Documents/projects/DO_project/results/flat/coloc/all_coloc_greaterorover75",quote = F,row.names = F)



## N.B. GPR133 and PTRF and ZNF609 are in nominal_hum as "Adgrd1" and "Cavin1" and "ZFP609", respectively
###
all_coloc_greaterorover75[which(all_coloc_greaterorover75$gene == "GPR133"),"gene"] = "ADGRD1"
all_coloc_greaterorover75[which(all_coloc_greaterorover75$gene == "PTRF"),"gene"] = "CAVIN1"
all_coloc_greaterorover75[which(tolower(all_coloc_greaterorover75$gene) == "ZNF609"),1] = "ZFP609"

###

#annotate kda tables with coloc data
kda_all = read.csv("./results/flat/key_driver_analysis_sexcombined_sft4.csv")
kda_F = read.csv("./results/flat/key_driver_analysis_FEMALES_sft4.csv")
kda_M = read.csv("./results/flat/key_driver_analysis_MALES_sft5.csv")

#
kda_all$eBMD = kda_F$eBMD = kda_M$eBMD = NA
kda_all$FNBMD = kda_F$FNBMD = kda_M$FNBMD = NA
kda_all$LSBMD = kda_F$LSBMD = kda_M$LSBMD = NA
kda_all$PPH4 = kda_F$PPH4 = kda_M$PPH4 = NA



#kda_all
for(i in 1:nrow(kda_all)){
  if(tolower(kda_all$gene[i]) %in% tolower(all_coloc_greaterorover75$gene) & kda_all$hyper[i]<=0.05){
    x = subset(all_coloc_greaterorover75, gene == toupper(kda_all$gene[i]))
    
    if("BMD" %in% x$pheno){kda_all$eBMD[i] = 1}
    if("FNBMD" %in% x$pheno){kda_all$FNBMD[i] = 1}
    if("LSBMD" %in% x$pheno){kda_all$LSBMD[i] = 1}
    kda_all$PPH4[i] = max(x$H4)
    
  }
}
########
#kda_F
#
for(i in 1:nrow(kda_F)){
  if(tolower(kda_F$gene[i]) %in% tolower(all_coloc_greaterorover75$gene) & kda_F$hyper[i]<=0.05){
    x = subset(all_coloc_greaterorover75, gene == toupper(kda_F$gene[i]))
    
    if("BMD" %in% x$pheno){kda_F$eBMD[i] = 1}
    if("FNBMD" %in% x$pheno){kda_F$FNBMD[i] = 1}
    if("LSBMD" %in% x$pheno){kda_F$LSBMD[i] = 1}
    kda_F$PPH4[i] = max(x$H4)
    
  }
}

########
#kda_M
#
for(i in 1:nrow(kda_M)){
  if(tolower(kda_M$gene[i]) %in% tolower(all_coloc_greaterorover75$gene) & kda_M$hyper[i]<=0.05){
    x = subset(all_coloc_greaterorover75, gene == toupper(kda_M$gene[i]))
    
    if("BMD" %in% x$pheno){kda_M$eBMD[i] = 1}
    if("FNBMD" %in% x$pheno){kda_M$FNBMD[i] = 1}
    if("LSBMD" %in% x$pheno){kda_M$LSBMD[i] = 1}
    kda_M$PPH4[i] = max(x$H4)
    
  }
}





#####ADD IMPC
kda_all$impc = kda_F$impc = kda_M$impc = NA


impc = read.csv("./results/flat/IMPC_BMD_Results.csv")
impc_gene = impc[which(impc$genotype_p_value<=0.05),"Gene"]

kda_all[which(tolower(kda_all$gene) %in% impc_gene),"impc"] = 1
kda_F[which(tolower(kda_F$gene) %in% impc_gene),"impc"] = 1
kda_M[which(tolower(kda_M$gene) %in% impc_gene),"impc"] = 1

#write.csv(kda_all, file = "~/Desktop/kda_combined.csv", quote = F, row.names = F)
#write.csv(kda_F, file = "~/Desktop/kda_F.csv", quote = F, row.names = F)

#write.csv(kda_M, file = "~/Desktop/kda_M.csv", quote = F, row.names = F)


