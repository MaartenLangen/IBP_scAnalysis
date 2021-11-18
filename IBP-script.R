library("BPSC")
library("SingleCellExperiment")
library("scater")

# Load the data (first 2500 genes, first 2500 cells) ---------------------------


# X..MatrixMarket references gene in excel file, matrix is replicate, coordinate
# is count
gene_by_cell_count_matrix <-
  read.csv("GSE126954_gene_by_cell_count_matrix.txt", sep = "")

# Only look at the first 2500 genes and 2500 cells.
gene_by_cell_count_matrix <- gene_by_cell_count_matrix[,1:3]
gene_by_cell_count_matrix <- subset(gene_by_cell_count_matrix, X..MatrixMarket <= 2500)
gene_by_cell_count_matrix <- subset(gene_by_cell_count_matrix, matrix <= 2500)

# generate matrix where row=gene, column=replicate
data_matrix = matrix(0,nrow=2500,ncol=2500)

for (i in 2:nrow(gene_by_cell_count_matrix)){
  geneID = gene_by_cell_count_matrix$X..MatrixMarket[i]
  repl = gene_by_cell_count_matrix$matrix[i]
  count = gene_by_cell_count_matrix$coordinate[i]
  data_matrix[geneID,repl]=count
}

# load gene and cell annotations
gene_annotation <- read.csv("GSE126954_gene_annotation.csv")[1:2500,2:3]
cell_annotation <- read.csv("GSE126954_cell_annotation.csv")[1:2500,-1]

row.names(data_matrix) <- gene_annotation[1:2500,1]
names(data_matrix) <- cell_annotation[1:2500,1]

# Cleaning the dataset ---------------------------------------------------------

## Remove genes that are not expressed in any cell -----------------------------
keep_feature <- rowSums(data_matrix > 0) > 0
data_matrix <- data_matrix[keep_feature,]

## Quality control: ------------------------------------------------------------
counts_per_cell = colSums(data_matrix)
hist(counts_per_cell, breaks = 100)

unique_genes_per_cell <- c(1:2500)

for (i in 1:ncol(data_matrix)){
  column = data_matrix[,i]
  unique_genes_per_cell[i] <- length(column[column>0])
}

hist(unique_genes_per_cell,breaks=100)

# TODO: make thresholds based on histograms

# TODO: normalize data

mat.res=estimateBPMatrix(data_matrix,para.num=4,fout=NULL,estIntPar=FALSE,useParallel=TRUE)
