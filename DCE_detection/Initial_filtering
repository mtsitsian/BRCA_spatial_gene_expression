####Initial filtering and setup of your count matrix

library(IOBR) #for filtering the duplicate isoforms in the count matrix ensembl ids
library(tidyverse) # For any data wrangling
library(edgeR) # Dependency for MDSeq package
library(MDSeq) # For count normalization
library(magrittr) # For pipeline execution
library(dendextend) # Check for unique values

#load the counts file from the TCGA summarized experiment, this is the count matrix you will need (gene names (preferrably ensemblids) in the rownames, subtype names as the column names)
counts <- read.csv("TCGA_raw_counts_GC.csv", header = T)

#in case you want to rename your columns accoriding to subtype name
names(counts)[startsWith(names(counts), 'CIN')] <- 'CIN'
names(counts)[startsWith(names(counts), 'MSI')] <- 'MSI'
names(counts)[startsWith(names(counts), 'GS')] <- 'GS'
names(counts)[startsWith(names(counts), 'EBV')] <- 'EBV'
names(counts)[startsWith(names(counts), 'Normal')] <- 'NP'

#table(colnames(counts)) check the number of cases, in this analysis the molecular subtypes
counts <- cbind("ensid" = row.names(counts), `row.names<-`(counts, NULL))
#make the first column the rownames of the dataset
#colnames(counts)[1] <- "ensid"
#counts <- counts %>% remove_rownames() %>% column_to_rownames( var = "ensid")

#remove any duplicate gene isoforms by converting the counts to the mean of the times the gene appears in the matrix
counts$ensid <- gsub("\\..*","",counts$ensid)
counts <- remove_duplicate_genes(eset = counts ,column_of_symbol = "ensid", method = "mean")
summary(duplicated(rownames(counts))) #check if there are duplicates, must return FALSE


#load the coordinates data and make it a dataframe
coordinates <- rtracklayer::import('Homo_sapiens.GRCh38.109.gtf') #https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/
coordinates=as.data.frame(coordinates)

#select the genes from the initial annotation data
coordinates <- coordinates[coordinates$type == "gene", ]

#change the "seqnames" column to "chromosome" and put "chr" in front of the numbers
coordinates$seqnames <- paste("chr", coordinates$seqnames, sep = '')
coordinates <- coordinates %>% rename("chromosome" = seqnames)

#make the gene_version column
coordinates$gene_version = paste(coordinates$gene_id, coordinates$gene_version, sep=".")

#select the columns you want to keep
coordinates <- coordinates %>% select(chromosome, start, end, strand, gene_id, gene_version, gene_name, gene_biotype)

#put the column of gene_id as rownames so you can merge them with the count matrix
coordinates <- coordinates %>% remove_rownames() %>% column_to_rownames(var = "gene_id")

#filter the coordinates to the preferred biotypes
coordinates <- coordinates[coordinates$gene_biotype %in% c("protein_coding",
                                                           "unprocessed_pseudogene",
                                                           "processed_pseudogene",
                                                           "lncRNA",
                                                           "transcribed_unprocessed_pseudogene",
                                                           "misc_RNA",
                                                           "transcribed_processed_pseudogene",
                                                           "transcribed_unitary_pseudogene",
                                                           "pseudogene"),]

table(coordinates$gene_biotype) #must include all the biotypes above

#filter the coordinates to the preferred chromosomes
coordinates <- coordinates[coordinates$chromosome %in% c("chr1","chr2","chr3","chr4",
                                                         "chr5","chr6", "chr7", "chr8",
                                                         "chr9","chr10","chr11","chr12",
                                                         "chr13","chr14","chr15","chr16",
                                                         "chr17","chr18","chr19","chr20",
                                                         "chr21","chr22","chrX","chrY"),]

unique(coordinates$chromosome) #check if you have kept the right chromosomes

#check if the ensids are all unique in the counts and coordinates data (must return TRUE)
all_unique(rownames(coordinates))
all_unique(rownames(counts))

#keep the genes that have 10 counts or more across all samples - low count filtering
keep <- rowSums(counts) >= 10
counts <- counts[keep,] 
all_unique(rownames(counts)) #56899 genes, all unique!, must return TRUE
any(is.na(rownames(counts))) #must return FALSE

#keep the genes that are present both in coordinates and the counts dataframe
coor_rows <- rownames(coordinates)
counts <- counts[coor_rows,] #55468 genes
counts <- na.omit(counts) #remove any NA rows, they will be a problem for the analysis
all_unique(rownames(counts)) #must return TRUE



#keep the coordinates of the genes after the low count filtering
nz_expression <- rownames(counts)
coordinates <- coordinates[nz_expression,]
any(is.na(rownames(coordinates))) #check if there are any NA values, must return FALSE
all.equal(rownames(coordinates), rownames(counts)) #must return TRUE, up to this point with filtering, genes are 50888, all unique!


#make the coordinates dataframe a list for further analysis 
coordinates <- list(coordinates = coordinates)
raw_counts_filtered <- as.matrix(counts)

#save the filtered counts in case you want to run differential expression analysis
#write.csv(raw_counts_filtered, "raw_counts_filtered_GC_ensids_latest.csv")

data_normalized <- normalize.counts(counts = raw_counts_filtered, method = "RLE", verbose = T)

# Normalization for gene length
lengths <- coordinates$coordinates$end - coordinates$coordinates$start + 1
lengths %<>% `/`(1000)
for (i in 1:nrow(data_normalized)){
  print(i)
  data_normalized[i,] <- data_normalized[i,]/lengths[i]
}

data_normalized <- as.data.frame(data_normalized)

names(data_normalized)[startsWith(names(data_normalized), 'CIN')] <- 'CIN'
names(data_normalized)[startsWith(names(data_normalized), 'MSI')] <- 'MSI'
names(data_normalized)[startsWith(names(data_normalized), 'GS')] <- 'GS'
names(data_normalized)[startsWith(names(data_normalized), 'EBV')] <- 'EBV'
names(data_normalized)[startsWith(names(data_normalized), 'NP')] <- 'NP'
table(colnames(data_normalized))
