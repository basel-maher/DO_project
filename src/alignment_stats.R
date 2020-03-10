set.seed(8675309)
options(digits = 10)

align_files = list.files(path = "./results/flat/RNA-seq/sums/", full.names = TRUE, recursive = F)

aln_tab = data.frame(matrix(ncol = 5), stringsAsFactors = F)
colnames(aln_tab) = c("run","sample","reads_processed","reads_w_atleast_1_aln","reads_fail2aln")
read_tab = data.frame(matrix(ncol=1,nrow=length(align_files)))

x  =read.delim(align_files[1],header = FALSE,stringsAsFactors = FALSE)
align_rate = strsplit(strsplit(x[15,],split = " ")[[1]][1],"%")[[1]][1]

reads = strsplit(x[1,],split = " ")[[1]][1]

for(f in 1:length(align_files)){
  
  x  =read.delim(align_files[f],header = FALSE,stringsAsFactors = FALSE)
  reads = strsplit(x[1,],split = " ")[[1]][1]
  print(reads)
  read_tab[f,1]=as.numeric(reads)
  rownames(read_tab)[f] = align_files[f]
}

samps = rownames(read_tab)
samp_names = c(rep(NA,times=length(samps)))

z = sapply(X = as.character(samps[c(1:48,97:192,289:384,433:816)]),function(x) unlist(strsplit(unlist(strsplit(x,split="_"))[5],split="-"))[1])
z = unname(z)
samp_names[c(1:48,97:192,289:384,433:816)] = z


z = sapply(X = as.character(samps[c(49:96,193:288)]),function(x) unlist(strsplit(unlist(strsplit(x,split="_"))[4],split="-"))[1])
z = unname(z)
samp_names[c(49:96,193:288)] = z

z = sapply(X = as.character(samps[c(385:432)]),function(x) unlist(strsplit(unlist(strsplit(x,split="_"))[3],split="-"))[1])
z = unname(z)
samp_names[c(385:432)] = z

read_tab$samp = samp_names

sort(aggregate(matrix.ncol...1..nrow...length.align_files..~samp,data=read_tab,FUN=sum)[,2])
#min = ~24.4 mil 
#max = ~55 mil 
#mean = ~39.2 mil
#median=~38.98mil

aln_tab = data.frame(matrix(ncol = 1,nrow=length(align_files)), stringsAsFactors = F)
for(f in 1:length(align_files)){
  x  =read.delim(align_files[f],header = FALSE,stringsAsFactors = FALSE)
  
  align_rate = strsplit(strsplit(x[15,],split = " ")[[1]][1],"%")[[1]][1]
  print(align_rate)
  aln_tab[f,1]=as.numeric(align_rate)
}

#overall alignment rates
#min=74.59
#max=96.04
#median=94.93
#mean=93.95

