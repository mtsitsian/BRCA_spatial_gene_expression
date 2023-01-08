######################DATA MINING FROM TCGA##########################

library(BiocManager)
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)

#Specify the parameters of the data you want to mine from TCGA
CancerProject <- "TCGA-BRCA"
dir <- "C:/Users/geots/OneDrive/Desktop/CG^2/TCGA-BRCA"

query <- GDCquery(project = CancerProject,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

samplesDown <- getResults(query,cols=c("cases"))

#download Primary tumor samples
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")

#download Normal tumor samples (Normal-like)
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")

queryDown <- GDCquery(project = CancerProject, 
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts", 
                      barcode = c(dataSmTP,dataSmNT))

GDCdownload(query = queryDown,directory = dir)                    

dataPrep <- GDCprepare(query = queryDown, save = TRUE,directory = dir, save.filename = "htseq_counts.rda", summarizedExperiment = TRUE)

####################### MINE GENE EXPRESSION DATA WITH THE BARCODES FROM TCGA##########################


listsamples <- strsplit(dataSmTP, ",")

# Query platform with a list of barcode 
query <- GDCquery(project = "TCGA-BRCA", 
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  barcode = listsamples)


# Download a list of barcodes 
GDCdownload(query)

# Prepare expression matrix with geneID in the rows and samples (barcode) in the columns
# rsem.genes.results as values
BRCARnaseqSE <- GDCprepare(query)

BRCAMatrix <- assay(BRCARnaseqSE) # or BRCAMatrix <- assay(BRCARnaseqSE,"raw_count")
