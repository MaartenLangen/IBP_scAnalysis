# Install the packages
# install.packages("devtools")
library("devtools")
# install_github("nghiavtr/BPSC")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("biomaRt")

# install.packages("remotes")
# remotes::install_github("davismcc/scater")

# Load the packages
library("BPSC")
library("SingleCellExperiment")
library("scater")
library(biomaRt)
library("Matrix")

# Load the data

# X..MatrixMarket references gene in excel file, matrix is replicate, coordinate
# is count
gene_by_cell_count_matrix <-
  read.csv("GSE126954_gene_by_cell_count_matrix.txt", sep = "")[,1:3]

# load gene and cell annotations
gene_annotation <- read.csv("GSE126954_gene_annotation.csv")[,2:3]
cell_annotation <- read.csv("GSE126954_cell_annotation.csv")[,-1]

cells_bin = cell_annotation[cell_annotation$raw.embryo.time.bin=="330-390",]

# generate matrix where row=gene, column=replicate
data_matrix = sparseMatrix(i=gene_by_cell_count_matrix$X..MatrixMarket,
                          j=gene_by_cell_count_matrix$matrix,
                          x=gene_by_cell_count_matrix$coordinate
                           )

rownames(data_matrix)=gene_annotation[,1]
colnames(data_matrix)=cell_annotation[,1]

# select a small part of data
data_matrix = data_matrix[,as.numeric(unlist(rownames(cells_bin)))][,1:1000]

# create SCE object
sce <- SingleCellExperiment(
  assays = list(counts = data_matrix),
  rowData = data.frame(gene_names = rownames(data_matrix)),
  colData = data.frame(cell_names = colnames(data_matrix))
)

# Cleaning the dataset ---------------------------------------------------------

## Remove genes that are not expressed in any cell -----------------------------
keep_feature <- rowSums(counts(sce) > 0) > 0
sce_filtered <- SingleCellExperiment(
  assays = list(counts = counts(sce)[keep_feature, ]),
  rowData = data.frame(gene_names = rowData(sce)[keep_feature,]),
  colData = data.frame(cell_names = colData(sce))
)

## Quality control: ------------------------------------------------------------
sce_filtered <- addPerCellQC(sce_filtered)
sce_filtered <- addPerFeatureQC(sce_filtered)

hist(sce_filtered$total, breaks = 100) #hist of counts per cell

unique_genes_per_cell <- colSums(counts(sce)>0)
hist(unique_genes_per_cell,breaks=100)

# TODO: make thresholds based on histograms

### Normalize data

# Use biomRT package to get access to WormBase ParaSite BioMart
# Establish a connection to the WormBase ParaSite BioMart
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# List available filters that can be used in a query and attributes you can retrieve
listFilters(mart)[1:100,]
listAttributes(mart)[1:100,]

# Retrieve the start and end positions of the genes, as well as the GO terms
gene_info <- getBM(filters="wbps_gene_id", 
  attributes=c("wbps_gene_id", "start_position","end_position","go_accession","go_name_1006"), 
  values=row.names(data_matrix), 
  mart=mart,
  uniqueRows=FALSE)
dim(gene_info) 
head(gene_info)
# biomaRt doesn't return NA if it can find it, so the size could not match to the data_matrix
# It also doesn't return results in the same order.

# Calculate the gene length
gene_info$length <- (gene_info$end_position - gene_info$start_position)
gene_info <- gene_info[-c(2:3)]
colnames(gene_info)[3] <- "go_name"
head(gene_info)




###iDEA
#calculate variance
#Here the first results is used!!
pvalue <- results[,2] #### the pvalue column
zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
fc <- results[,1] ## the fold change column
se_beta <- abs(fc/zscore) ## to approximate the standard error of beta
var = se_beta^2  ### square 
summary = data.frame(fc = fc,variance = var)# Summary is a matrix of fold change and variance of each gene
