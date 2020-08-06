
library(PhenStat)

## creates a with weight and without weight results table
## set x to any table with a getDataFunction, Gene, PhenCenter.and Pipeline column
DEXATable <- read.csv("IMPC_DEXA_Table.csv")
x <- DEXATable
## x <- SampleNtable # go back later to check if more samples were added
## x <- Not2Genotypes # go back later to check if appropriate control / test strains were added

###############################################################

for(i in 1:nrow(x)) {
  g <- eval(parse(text=as.character(x$getDataFunction[i])))
  g$Weight <-as.numeric(as.character(g$Weight))
  g$Genotype <- g$biological_sample_group
  g$biological_sample_group <- NULL
  cCheck <- subset(g, Genotype == "control", select= c(Sex, Zygosity))
  xCheck <- subset(g, Genotype == "experimental", select= c(Sex, Zygosity))
  het <- sum(xCheck$Zygosity == "heterozygote")
  homo <- sum(xCheck$Zygosity == "homozygote")
  hemi <- sum(xCheck$Zygosity == "hemizygote")
  zygo <- paste(het, "-hets/", homo, "-homos/", hemi, "-hemis",sep ='')
  gene <- as.character(x$Gene[i])
  phenCenter <- as.character(x$PhenCenter[i])
  pipeline <- as.character(x$Pipeline[i])
  getDataFunction <- as.character(x$getDataFunction[i])
  
  if(nlevels(unique(g$Genotype)) < 2){
    TempResults[1,] <- "Not 2 Genotypes"
    TempResults2[1,] <- "Not 2 Genotypes"
  } else if (sum(cCheck$Sex == "male") < 4
             || sum(cCheck$Sex =="female") < 4
             || sum(xCheck$Sex == "male") < 4
             || sum(xCheck$Sex == "female") < 4){
    TempResults[1,] <- "Sample N"
    TempResults2[1,] <- "Sample N"
  } else {
    test <- PhenList(dataset= g, testGenotype= "experimental", refGenotype="control")
    result <- testDataset(test, depVariable="Value", equation="withWeight")
    result2 <- testDataset(test, depVariable="Value", equation="withoutWeight")
    TempResults <- as.data.frame(t(result@analysisResults$model.output.summary))
    TempResults2 <- as.data.frame(t(result2@analysisResults$model.output.summary))
    TempResults$genotype_p_value <- result@analysisResults[["model.output.genotype.nulltest.pVal"]]
    TempResults2$genotype_p_value <- result2@analysisResults[["model.output.genotype.nulltest.pVal"]]
  }
  
  TempResults$X=NULL
  TempResults$Gene <- gene
  TempResults$PhenCenter <- phenCenter
  TempResults$Pipeline <- pipeline
  TempResults$getDataFunction <- getDataFunction
  TempResults$test_zygosity <- zygo
  TempResults2$X=NULL
  TempResults2$Gene <- gene
  TempResults2$PhenCenter <- phenCenter
  TempResults2$Pipeline <- pipeline
  TempResults2$getDataFunction <- getDataFunction
  TempResults2$test_zygosity <- zygo
  
  col_idx <- grep("test_zygosity", names(TempResults))
  TempResults <- TempResults[, c(col_idx, (1:ncol(TempResults))[-col_idx])] 
  col_idx <- grep("Gene", names(TempResults))
  TempResults <- TempResults[, c(col_idx, (1:ncol(TempResults))[-col_idx])]
  col_idx <- grep("test_zygosity", names(TempResults2))
  TempResults2 <- TempResults2[, c(col_idx, (1:ncol(TempResults2))[-col_idx])] 
  col_idx <- grep("Gene", names(TempResults2))
  TempResults2 <- TempResults2[, c(col_idx, (1:ncol(TempResults2))[-col_idx])]
  
  if(i == 1) {
    ResultsTable <- TempResults
    ResultsTable2 <- TempResults2
  } else {
    ResultsTable <-  rbind(ResultsTable, TempResults)
    ResultsTable2 <-  rbind(ResultsTable2, TempResults2)
  }
  if(i%%10 == 0) {
    write.csv(ResultsTable,"IMPC_BMD_Results")
    View(ResultsTable)
    write.csv(ResultsTable2,"IMPC_BMD_Results_noWeight")
    View(ResultsTable2)
  }
}
##################################################################################


## Table of "Sample N" error genes using subset function
SampleNtable <- subset(ResultsTable, genotype_estimate == "Sample N",
                       select=c(Gene, PhenCenter, Pipeline, getDataFunction))
write.csv(SampleNtable,"Sample_N_Errors")

## Table of "Not 2 Genotypes" error genes using subset function
Not2Genotypes <- subset(ResultsTable, genotype_estimate == "Not 2 Genotypes",
                        select=c(Gene, PhenCenter, Pipeline, getDataFunction))
write.csv(Not2Genotypes,"Not_2_Genotypes")

## create final resuls table filtered of error rows
FinalResults <- ResultsTable[-c(which(ResultsTable$genotype_estimate_SE == "Sample N")), ] 
FinalResults <- FinalResults[-c(which(FinalResults$genotype_estimate_SE == "Not 2 Genotypes")), ] 
FinalResult$Gene <- tolower(FinalResults$Gene)
FinalResults <- FinalResults[order(FinalResults$genotype_p_value),]

## not sexually dimorphic table
ResultsTable$genotype_estimate <-as.numeric(as.character(ResultsTable$genotype_estimate))
NotDimorphic <- subset(ResultsTable, genotype_estimate > -100 & (genotype_p_value * nrow(FinalResults)) < 0.05,
                       select=c(Gene, genotype_estimate, genotype_estimate_SE, genotype_p_value, sex_estimate,
                                sex_estimate_SE, sex_p_value, weight_estimate, weight_estimate_SE, weight_p_value,
                                intercept_estimate, intercept_estimate_SE, PhenCenter, Pipeline, getDataFunction))
write.csv(NotDimorphic,"Not_Dimorphic_Results.csv")

##subset sexual dimorphs
ResultsTable$sex_MvKO_p_value <- as.numeric(as.character(ResultsTable$sex_MvKO_p_value))
ResultsTable$sex_FvKO_p_value <- as.numeric(as.character(ResultsTable$sex_FvKO_p_value))
Dimorphs <- subset(ResultsTable, sex_MvKO_p_value >= 0,
                   select=c(Gene, sex_FvKO_p_value, sex_FvKO_estimate, sex_FvKO_SE, sex_MvKO_p_value, sex_MvKO_estimate,
                            sex_MvKO_SE, genotype_p_value, sex_estimate, sex_estimate_SE, sex_p_value, weight_estimate, 
                            weight_estimate_SE, weight_p_value, intercept_estimate, intercept_estimate_SE, PhenCenter,
                            Pipeline, getDataFunction))

## significant male dimorphs
MaleDimorphs <- subset(Dimorphs, sex_MvKO_p_value < sex_FvKO_p_value & (sex_MvKO_p_value * nrow(FinalResults)) < 0.05,
                       select=c(Gene, sex_MvKO_p_value, sex_MvKO_estimate, sex_MvKO_SE, sex_FvKO_p_value, genotype_p_value, sex_estimate, sex_estimate_SE, sex_p_value, weight_estimate, 
                                weight_estimate_SE, weight_p_value, intercept_estimate, intercept_estimate_SE, PhenCenter,
                                Pipeline, getDataFunction))
MaleDimorphs <- MaleDimorphs[order(MaleDimorphs$sex_MvKO_p_value),]
rownames(MaleDimorphs) <- 1:nrow(MaleDimorphs)
write.csv(MaleDimorphs,"Male_Dimorph_Results.csv")

## significant female dimorphs
FemaleDimorphs <- subset(Dimorphs, sex_FvKO_p_value < sex_MvKO_p_value & (sex_FvKO_p_value * nrow(FinalResults)) < 0.05,
                         select=c(Gene, sex_FvKO_p_value, sex_FvKO_estimate, sex_FvKO_SE, sex_MvKO_p_value, genotype_p_value, sex_estimate, sex_estimate_SE, sex_p_value, weight_estimate, 
                                  weight_estimate_SE, weight_p_value, intercept_estimate, intercept_estimate_SE, PhenCenter,
                                  Pipeline, getDataFunction))
FemaleDimorphs <- FemaleDimorphs[order(FemaleDimorphs$sex_FvKO_p_value),]
rownames(FemaleDimorphs) <- 1:nrow(FemaleDimorphs)
write.csv(FemaleDimorphs,"Female_Dimorph_Results.csv")
