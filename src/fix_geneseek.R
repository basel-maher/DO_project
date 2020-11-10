#Fix geneseek files
#library(qtl2convert) dont think i need this

# - There was a sample confusion. Re-genotyped samples have a .1 appended to their name, except 371.
# - In ./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport, change Sample ID 371 to 371.1

#skip 9 to skip the header lines in FinalReports
myfread_fr <- function(filename) data.table::fread(filename, data.table=FALSE,skip = 9)


###change 371 to 371.1 in /Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport.txt"
x = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport.txt")
x$`Sample ID`[which(x$`Sample ID`=="371")] = "371.1"

write.table(x = x, file = "./data/GIGAMUGA/FinalReport_files/20190228_FinalReport_FIXED.txt", quote = FALSE,row.names = FALSE,sep = "\t")
rm(x)
#MAKE SURE TO ADD HEADER BACK MANUALLY!!

#Header is the following, without hashes
# [Header]
# GSGT Version	2.0.2
# Processing Date	2/28/2019 11:49 AM
# Content		GigaMuga_11769261_A.bpm
# Num SNPs	143259
# Total SNPs	143446
# Num Samples	96
# Total Samples	96
# [Data]


# - Fix ./data/GIGAMUGA/merged/Merged_Sample_Map.txt manually. Change the second 371 entry to 371.1

# create merged_FinalReport.txt

# - Generate ./data/GIGAMUGA/merged/Merged_FinalReport.txt.
#These are GeneSeek FinalReport files
x = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20171110_FinalReport.txt")
x2 = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20180408_FinalReport.txt")
x3 = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20181213_FinalReport.txt")
x4 = myfread_fr("./data/GIGAMUGA/FinalReport_files/20190228_FinalReport_FIXED.txt")


merged = rbind(x,x2,x3,x4)

write.table(x = merged, file = "./data/GIGAMUGA/merged/Merged_FinalReport.txt", quote = FALSE,row.names = FALSE,sep = "\t")

rm(x,x2,x3,x4,merged)
###
#ADD THIS HEADER TO THE FILE MANUALLY!! without hashes

# [Header]
# GSGT Version	2.0.2
# Processing Date	11/14/2017 10:35 AM
# Content		GigaMuga_11769261_A.bpm
# Num SNPs	143259
# Total SNPs	143446
# Num Samples	744
# Total Samples	744
# [Data]

