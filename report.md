# Problem statement

Using the publicly available bulk RNA-Seq data from GTEX (v8), identify biological pathways/processes specific to the liver. As what constitutes a biological pathway/process is open for interpretation, please use your preferred method / definition / approach (or multiple). For bonus points, identify putative transcriptional regulators of these pathways/processes. To ensure that you do not encounter memory issues please restrict your analyses to protein-coding genes and the following tissues: heart, kidney, liver, lung, muscle, pancreas, spleen, stomach, pituitary gland, and thyroid.

# Obtaining data
Simple download of files from GTEx website

# Preprocessing data
## Find protein coding genes
Using gencode file
Distribution of lncRNA is different from other genes

## Subset to tissues of interest

## Subset to samples of adequate quality

## Identify and correct metadata errors

## Identify RNA contamination
(gene from other samples leaking) [ref]

# Identifying pathways specific to the liver
## Definition of specificity
"Specific" pathway could be interpreted as the task of identifying:
- differential expression, i.e., a differential increase or decrease in activity compared to other tissues
- more strictly, presence or absence compared to other tissues (similar to the task of identifying "markers")

## Definition of pathway
Here we will simply define pathways based on pre-defined gene sets present in online databases. 

## Analysis strategy

### Approach 0: expression
- For each tissue, identify genes expressed or not
- Take the liver-specific ones
- Perform GSEA analysis

### Approach 1: single sample
- For each tissue, identify co-expressed gene modules using WGCNA
- For each tissue, perform GSEA. 
- Take the liver-specific pathways

### Approach 2: single sample
- For each tissue, perform ssGSEA

### Approach 3
- Find DE genes from testing: 
    - all vs liver
    - pairwise vs liver

### Approach 4
Bayesian approach from GTEx pilot


# Identify transcriptional regulators


# References
GTEx pilot study
yarn normalization
GSEA
ORA
Progeny