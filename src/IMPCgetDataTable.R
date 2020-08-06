if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("PhenStat")
BiocManager::install("IMPCdata")
library(IMPCdata)
library(PhenStat)

## Print all phenCenters with their respective pipelines
phencens <- unlist(getPhenCenters())
for(center in phencens){  
  print(center)
  print(getPipelines(center))
}
  
  ## Example of how to get all genes with DEXA data from a known phenCenter, "TCP"
  getPipelines("TCP")
  
  ##Next check phenCenter, pipeline to confirm procedure "IMPC_DXA_001" is returned
  getProcedures("TCP","TCP_001")
  
  ## using get table to return table with functions to call data from DEXA procedure at TCP
    ## first argument is name of returned csv file
    ## 2nd arg is phenCenter, 3rd is pipeline, 4th is procedure, 5th is parameter)
  getIMPCTable("IMPC_TCP_DEXA", "TCP", "TCP_001", "IMPC_DXA_001", "IMPC_DXA_004_001" )
  
  
  ## Defining function to call above functions sequentially and merge tables into a single table.csv
  GetIMPCDEXATable <- function() {
    getIMPCTable("IMPC_TCP_DEXA", "TCP", "TCP_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    DEXATable <- read.csv(file="IMPC_TCP_DEXA_1.csv", header=TRUE, sep=",")
    getIMPCTable("IMPC_JAX_DEXA", "JAX", "JAX_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_JAX_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable)
    getIMPCTable("IMPC_WTSI_DEXA", "WTSI", "MGP_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_WTSI_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable)
    getIMPCTable("IMPC_BCM_DEXA", "BCM", "BCM_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_BCM_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable) 
    getIMPCTable("IMPC_HMGU_DEXA", "HMGU", "HMGU_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_HMGU_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable)  
    getIMPCTable("IMPC_ICS_DEXA", "ICS", "ICS_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_ICS_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable)  
    getIMPCTable("IMPC_KMPC_DEXA", "KMPC", "IMPC_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_KMPC_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable)  
    getIMPCTable("IMPC_MARC_DEXA", "MARC", "IMPC_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_MARC_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable) 
    getIMPCTable("IMPC_HRWL_DEXA", "MRC Harwell", "HRWL_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_HRWL_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable) 
    getIMPCTable("IMPC_UCD_DEXA", "UC Davis", "UCD_001", "IMPC_DXA_001", "IMPC_DXA_004_001")
    TempTable <- read.csv(file="IMPC_UCD_DEXA_1.csv", header=TRUE, sep=",")
    DEXATable <-  rbind(DEXATable, TempTable) 
    write.csv(DEXATable,"IMPC_DEXA_Table.csv")  }
  ##Returns individual DEXA Tables for each PhenCenter and a total merged table: IMPC_DEXA_Table.csv
  GetIMPCDEXATable()
  
  ## DEXATable is used with IMPCPhenStats script to get data and run through phenstats