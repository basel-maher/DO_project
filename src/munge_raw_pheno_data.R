##This script imports raw data from collaborators, outputs pheno table##
##raw data in /data/pheno_data/

##############################################
##Bending data
#read raw data
bending_data = read.csv("./data/pheno_data/DO mice 4-point bending.csv", stringsAsFactors = FALSE)


#remove "DO-" from identifier
bending_data$Specimen = apply(bending_data,1,function(x) sub("DO-","",x[1]))

#edit colnames
colnames(bending_data) = c("specimen","femur_length","stiffness","yield_stiffness","max_load","yield_load",
                           "frax_load","PYD","disp_at_yield","disp_at_frax","disp_at_max_load",
                           "total_work","work_to_yield","work_post_yield")
#add "bending_" to colnames
colnames(bending_data) <- paste("bending", colnames(bending_data), sep = "_")


##############################################
##Histo data
#read raw data
histo_data = read.csv("./data/pheno_data/histo_master.csv", stringsAsFactors = FALSE)

#remove "GROUP" column
histo_data = histo_data[,-2]

#remove "FV" and 'FV0" from identifier in SPEC
histo_data$SPEC.. = apply(histo_data,1,function(x) sub("FV0","",x[1]))
histo_data$SPEC.. = apply(histo_data,1,function(x) sub("FV","",x[1]))

#rename 1st col
colnames(histo_data)[1] = c("specimen")

#add "histo_" to columns
colnames(histo_data) <- paste("histo", colnames(histo_data), sep = "_")

##############################################
##uCT data
#read raw data
microCT_data = read.csv("./data/pheno_data/uCT_master.csv",stringsAsFactors = FALSE)

#change colnames
colnames(microCT_data) = c("specimen","sex","batch","femur_length","BV.TV","BMD","BS.BV","conn_density","SMI","Tb.N",
                           "Tb.Th","Tb.Sp","Ct.Ar","Ma.Ar","Tt.Ar","Ct.Ar/Tt.Ar","Ct.Th","Ct.TMD","Ct.porosity", "pMOI","Imax","Imin")

#add uCT_ to colnames"
colnames(microCT_data) <- paste("uCT", colnames(microCT_data), sep = "_")

####
#EVEN MORE RAW DATA. IGNORE FOR NOW#
# microCT_data_trab = read.csv("uCT_trab_1-578.csv", stringsAsFactors = FALSE)
# colnames(microCT_data_trab)[3]="uCT_batch"
# colnames(microCT_data_trab) <- paste("uCT_trab", colnames(microCT_data_trab), sep = "_")
# 
# microCT_data_cort_out = read.csv("uCT_cort_out_1-578.csv", stringsAsFactors = FALSE)
# colnames(microCT_data_cort_out)[3]="uCT_batch"
# microCT_data_cort_out = (na.omit(microCT_data_cort_out))
# colnames(microCT_data_cort_out) <- paste("uCT_cort_out", colnames(microCT_data_cort_out), sep = "_")
# 
# microCT_data_cort_donut = read.csv("uCT_cort_donut_1-578.csv", stringsAsFactors = FALSE)
# colnames(microCT_data_cort_donut)[3]="uCT_batch"
# colnames(microCT_data_cort_donut) <- paste("uCT_cort_donut", colnames(microCT_data_cort_donut), sep = "_")
# 
#####

##############################################
#MAT data
#read raw data
MAT_data = read.csv("./data/pheno_data/MAT_data.csv", stringsAsFactors = FALSE)

#remove column 6, contained notes
MAT_data = MAT_data[,-6]

#set the column with * in it to NA - specimen cut short
MAT_data[which(MAT_data$VOI4=="*"),"VOI4"] = NA

#remove "FV" and "FV0" from sample names
MAT_data$sample.. = apply(MAT_data,1,function(x) sub("FV0","",x[1]))
MAT_data$sample.. = apply(MAT_data,1,function(x) sub("FV","",x[1]))

#rename cols
colnames(MAT_data) = c("specimen","VOL1","VOL2","VOL3","VOL4")

#add MAT_ to colnames
colnames(MAT_data) <- paste("MAT", colnames(MAT_data), sep = "_")
#####
##############################################

##read in harvest records
harvest_data = read.csv("./data/pheno_data/Harvest Records.csv", stringsAsFactors = FALSE)

#subset to most recently completed rows
harvest_data = harvest_data[1:698,]

#create full_pheno_table by merging data 
full_pheno_table = merge(harvest_data, bending_data, by.x="Mouse.ID", by.y = "bending_specimen", all = T,sort = FALSE)

full_pheno_table = merge(full_pheno_table,histo_data,by.x = "Mouse.ID",by.y = "histo_specimen",all = T,sort = FALSE)

full_pheno_table = merge(full_pheno_table, microCT_data, by.x="Mouse.ID", by.y = "uCT_specimen", all = T,sort = FALSE)

full_pheno_table = merge(full_pheno_table, MAT_data, by.x="Mouse.ID", by.y = "MAT_specimen", all = T,sort = FALSE)

#READ IN OTHER uCT DATA. IGNORE#
#full_pheno_table = merge(full_pheno_table, microCT_data_cort_donut, by.x="Mouse.ID", by.y = "uCT_cort_donut_Sample.Name", all = T,sort = FALSE)
#full_pheno_table = merge(full_pheno_table, microCT_data_cort_out, by.x="Mouse.ID", by.y = "uCT_cort_out_Sample.Name", all = T,sort = FALSE)
#full_pheno_table = merge(full_pheno_table, microCT_data_trab, by.x="Mouse.ID", by.y = "uCT_trab_Sample.Name", all = T,sort = FALSE)
####

##############################################

#order rows by Mouse.ID
full_pheno_table = full_pheno_table[order(as.numeric(full_pheno_table$Mouse.ID)),]

#remove mice used for breeding
#20 mice
full_pheno_table = full_pheno_table[-which(full_pheno_table$Comments == "used for breeding"),]

#replace commas with semicolons
full_pheno_table <- sapply(full_pheno_table, gsub, pattern = ",", replacement= ";")
full_pheno_table = as.data.frame(full_pheno_table,stringsAsFactors = FALSE)


#remove unnecessary cols
#removed bending_femur_length, uCT_sex, uCT_batch, uCT_femur_length

full_pheno_table = full_pheno_table[,-c(23,54:56)] 

#rename some cols
colnames(full_pheno_table)[2:22] = c("sac_date","coat_color","sac_time","sex","DOB","age_at_sac_days","body_weight","body_length",
                                     "glucose","RFP","GFP","BFP","FFP","adiposity","comments","soleus_weight","gastroc_weight","FL",
                                     "ML","AP","DO_generation")

#recalculate adiposity. It is the sum of RFP, GFP, and FFP in mg, divided by BW in mg, times 100 (ratio)
full_pheno_table$adiposity = ((as.numeric(full_pheno_table$RFP) + as.numeric(full_pheno_table$GFP) + as.numeric(full_pheno_table$FFP))
                              / (as.numeric(full_pheno_table$body_weight)*1000)) * 100



# because MAT=0 and MAT=NA are different, and because the distribution is zero-inflated, create two new "phenotypes" for each MAT pheno
#Namely, one which is only the non-zero values, and one that is binary (bin) for the presence or absence of MAT (nonzero) (this will be mapped as a binary trait)


full_pheno_table$MAT_VOL1_bin = full_pheno_table$MAT_VOL1_nonzero = full_pheno_table$MAT_VOL1 #create bin and nonzero phenos
full_pheno_table[which(full_pheno_table$MAT_VOL1 >0),"MAT_VOL1_bin"] = 1 #if greater than zerp, bin pheno = 1
full_pheno_table[which(full_pheno_table$MAT_VOL1 ==0),"MAT_VOL1_nonzero"] = NA # if zero, nonzero is NA

full_pheno_table$MAT_VOL2_bin = full_pheno_table$MAT_VOL2_nonzero = full_pheno_table$MAT_VOL2 #create bin and nonzero phenos
full_pheno_table[which(full_pheno_table$MAT_VOL2 >0),"MAT_VOL2_bin"] = 1 #if greater than zerp, bin pheno = 1
full_pheno_table[which(full_pheno_table$MAT_VOL2 ==0),"MAT_VOL2_nonzero"] = NA # if zero, nonzero is NA

full_pheno_table$MAT_VOL3_bin = full_pheno_table$MAT_VOL3_nonzero = full_pheno_table$MAT_VOL3 #create bin and nonzero phenos
full_pheno_table[which(full_pheno_table$MAT_VOL3 >0),"MAT_VOL3_bin"] = 1 #if greater than zerp, bin pheno = 1
full_pheno_table[which(full_pheno_table$MAT_VOL3 ==0),"MAT_VOL3_nonzero"] = NA # if zero, nonzero is NA

full_pheno_table$MAT_VOL4_bin = full_pheno_table$MAT_VOL4_nonzero = full_pheno_table$MAT_VOL4 #create bin and nonzero phenos
full_pheno_table[which(full_pheno_table$MAT_VOL4 >0),"MAT_VOL4_bin"] = 1 #if greater than zerp, bin pheno = 1
full_pheno_table[which(full_pheno_table$MAT_VOL4 ==0),"MAT_VOL4_nonzero"] = NA # if zero, nonzero is NA




#convert numerics to numeric
for(i in c(1,7:15,17:21,23:82)){
  full_pheno_table[,i] = as.numeric(full_pheno_table[,i])
}
##############################################

write.table(full_pheno_table,"./results/flat/full_pheno_table.csv",quote = FALSE, sep = ",",row.names = FALSE )
save(full_pheno_table, file = "./results/Rdata/full_pheno_table.Rdata")
