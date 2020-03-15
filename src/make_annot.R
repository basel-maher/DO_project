#make sure there are no factors. factors mess everything up and are a terrible invention
options(stringsAsFactors = FALSE)

#make a unified annotation file from the gene abundance files (stringtie output)
f = list.files("./results/flat/RNA-seq/abund/")

z=data.frame()

for(i in f){
  path = paste0("./results/flat/RNA-seq/abund/",i)
  tab = read.delim(path)
  
  z= rbind(z,tab)
  print(i)
}


z = unique(z)
#save(z, file="~/Desktop/z.RData")

#z_trim = z[-which(z$TPM == 0),]
z_trim=z
#take largest transcript as main
z_trim2=as.data.frame(matrix(ncol=ncol(z_trim), nrow=100000))

#keep longest transcript
colnames(z_trim2) = colnames(z_trim)
z_trim2$len = NA
i=1
for(transc in unique(z_trim$Gene.ID)){
  x = subset(z_trim, z_trim$Gene.ID== transc)
  x$len = abs(as.numeric(x$End) - as.numeric(x$Start))
  x = x[which(x$len == max(x$len)),][1,]
  #x = x[,-6]
  z_trim2[i,] =  x
  i = i+1
  print(i)
}


z = z_trim2[-which(is.na(z_trim2$Reference)),]

chr = c(1:19, "X")

z = z[-which(z$Reference %in% chr == FALSE),]

##there are 89 duplicated gene names
z[which(duplicated(z$Gene.Name)),]

#for now, give them a unique ID

z[which(duplicated(z$Gene.Name)),"Gene.Name"] = paste0(z[which(duplicated(z$Gene.Name)),"Gene.Name"], "_isoform")


which(duplicated(z$Gene.Name))

z$Gene.Name[7748] = paste0(z$Gene.Name[7748],".2")
z$Gene.Name[11056] = paste0(z$Gene.Name[11056],".2")
z$Gene.Name[34526] = paste0(z$Gene.Name[34526],".2")
z$Gene.Name[34527] = paste0(z$Gene.Name[34527],".3")
z$Gene.Name[48778] = paste0(z$Gene.Name[48778],".2")


which(duplicated(z$Gene.Name))

#z = z[,c(1,2)]
z = z[,-z$len]
write.csv(z,file="./results/flat/annot_file.csv", row.names = F, quote = F)



