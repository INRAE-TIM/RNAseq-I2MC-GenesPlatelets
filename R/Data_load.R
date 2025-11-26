# Loading assets #################################################################

source("./R/install_packages.R")
source('./R/custom_functions.R')

## Loading data ####
### pData (Samples file) ####
MetaData <- read.xlsx("./Data/Compte plaquettaire RNA Seq.xlsx",
                   sep = "\t")
pData <- data.frame(group = gsub("\\.(\\d+)", ".REP\\1", MetaData$Diet_Treatment.Rep))
row.names(pData) <- MetaData$Diet_Treatment.Rep
pData$group <- gsub('.REP[0-9]', '', pData$group)
pData <- cbind.data.frame(pData, 
                          do.call(rbind,
                                  strsplit(pData$group, '_')))
colnames(pData) <- c('group', 'diet', 'treatment')
pData$Rep <- do.call(rbind, strsplit(row.names(pData), '\\.'))[,2]
pData$TP <- MetaData$TP


### WorkingSet (Counts file) ####
Counts <- read.table(
  "./Data/gene_expression.xls",
  sep = "\t",
  header = T
) 
colnames(Counts)[colnames(Counts) == 'gene_symbol'] <- 'gene_name'
Counts <- Counts[,colnames(Counts)!='st_gene_id']

# Fitting data to the correct format ###########################################

## Ordering pData to the same order that Counts ####
pData <- pData[colnames(Counts)[grepl("\\d+", colnames(Counts))], ]

## Evaluating missing TP with missMDA ####

no.miss <- missMDA::imputePCA(
  rbind.data.frame(Counts[, grepl('[0-9]', colnames(Counts))], pData$TP),
  ncp = 5
  )$completeObs

colnames(no.miss)

MetaData$TP[is.na(MetaData$TP)]<- no.miss[dim(no.miss)[1], "HFC_CT.4"]
MetaData <- MetaData[!MetaData$Diet_Treatment.Rep %in% c('HFC_CT.1', 'HFC_THR.1'),]
MetaData$Diet_Treatment.Rep <-  gsub('CT', 'CTL', gsub('THR', 'LPC', gsub(
  'ROMI', 'HPC', MetaData$Diet_Treatment.Rep
)))

Counts <- Counts[,!colnames(Counts) %in% c('HFC_CT.1', 'HFC_THR.1')]
colnames(Counts) <- gsub('CT', 'CTL', gsub('THR', 'LPC', gsub(
  'ROMI', 'HPC', colnames(Counts)
)))

pData <- pData[!row.names(pData) %in% c('HFC_CT.1', 'HFC_THR.1'),]
row.names(pData) <- gsub('CT', 'CTL', gsub('THR', 'LPC', gsub(
  'ROMI', 'HPC', row.names(pData)
)))

pData$group <- gsub('CT', 'CTL', gsub('THR', 'LPC', gsub(
  'ROMI', 'HPC', pData$group
)))

pData$treatment <- gsub('CT', 'CTL', gsub('THR', 'LPC', gsub(
  'ROMI', 'HPC', pData$treatment
)))

rm(no.miss)