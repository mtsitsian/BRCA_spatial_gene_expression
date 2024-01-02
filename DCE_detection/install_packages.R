#Install the packages for the pipeline and load them in your working environment
#Grab a coffee or go for a walk, this will take long...

#Install the Bioconductor supported packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("tidyverse", "magrittr", "TCGAbiolinks", "SummarizedExperiment", "parallel", "edgeR", "dendextend",
"ComplexHeatmap", "features", "feather", "rlist", "valr", "venn", "clusterProfiler", "DOSE", "enrichplot", "ChIPpeakAnno", "AnnotationDbi", "org.Hs.eg.db", "ggsci", "ggplot2", "plyr",
"plotgardener", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg38.knownGene", "ggpubr", "gginnards","pheatmap", "karyoploteR", "DESeq2", "preprocessCore","cqn"
))

#Install IOBR and MDSeq packages
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")  
#BiocManager::install("ComplexHeatmap")

install.packages("doParallel")
install.packages("remotes")
install.packages("shiny")
install.packages("shinythemes")


#install IOBR
# options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
# options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

depens<-c('tibble', 'survival', 'survminer', 'sva', 'limma', "DESeq2","devtools",
          'limSolve', 'GSVA', 'e1071', 'preprocessCore', 'ggplot2', "biomaRt",
          'ggpubr', "devtools", "tidyHeatmap", "caret", "glmnet", "ppcor", "timeROC","pracma")
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))
    BiocManager::install(depen,update = FALSE)
}
remotes::install_github("IOBR/IOBR", force = T)
   
library(devtools)
options(download.file.method = "wininet")
install_github("zjdaye/MDSeq")

library(BiocManager)
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(parallel)
library(IOBR)
library(tidyverse)
library(edgeR)
library(MDSeq)
library(magrittr)
library(dendextend)
library(ComplexHeatmap)
library(features)
library(feather)
library(rlist)
library(valr)
library(venn)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(ChIPpeakAnno)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggsci)
library(ggplot2)
library(plyr)
library(plotgardener)
library(AnnotationHub)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggpubr)
library(gginnards)
library(pheatmap)
library(karyoploteR)
library(doParallel)
library(shiny)
library(shinythemes)
