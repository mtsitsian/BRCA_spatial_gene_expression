# Spatial gene expression analysis in Breast Cancer
Repository for the analysis of spatial gene expression in molecular subtypes of Breast Cancer

This is the repository for the code and data produced during the analysis of spatial organization of gene expression, performed in Breast Cancer transciptome profiles.

* [Description of the data](#description) 
* [Description of the code](#code)

## Description of the data<a name="description"></a>
The association of neighbouring genes at the expression level and the subsequent modelling of this at the level of whole genome topology have been described previously. In the present study, the same type of analysis performed in [Systemic Lupus Erythematosus](https://github.com/mtsitsian/SLE_spatial_gene_expression) is applied, respectively this time to a new transcriptomics dataset involving breast cancer subtypes. After each chromosome was segmented into bins of equal length, the correlation between their expression values was calculated and followed by the definition of regions of statistically significant coexpression events by applying a permutation test between local minima. The Domains of Coordinated Expression (DCEs) obtained from the analysis were verified by permutation analysis, which involves randomly redistributing the expression values of the bins (x 1000), for each chromosome. In this way, a statistically robust definition of chromosomal regions within which gene co-expression is statistically higher than in either region is achieved.


## Description of the code<a name="code"></a>
The analysis has been performed in R. All the code used to define, 
detect and analyze DCEs can be found inside the [DCE_Analysis](https://github.com/mtsitsian/BRCA_spatial_gene_expression/tree/main/DCE_analysis) and the [Dce_Detection](https://github.com/mtsitsian/BRCA_spatial_gene_expression/tree/main/DCE_detection) folder.
