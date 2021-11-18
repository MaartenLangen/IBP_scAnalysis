# Install the packages
install.packages("devtools")
library("devtools")
install_github("nghiavtr/BPSC")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
BiocManager::install("biomaRt")

install.packages("remotes")
remotes::install_github("davismcc/scater")

# Load the packages
library("BPSC")
library("SingleCellExperiment")
library("scater")
library(biomaRt)

# Load the data (first 2500 genes, first 2500 cells) ---------------------------


# X..MatrixMarket references gene in excel file, matrix is replicate, coordinate
# is count
gene_by_cell_count_matrix <-
  read.csv("C:/Users/oysdfx/Desktop/IBP/GSE126954_gene_by_cell_count_matrix.txt.gz", sep="")

# Only look at the first 2500 genes and 2500 cells.
gene_by_cell_count_matrix <- gene_by_cell_count_matrix[,1:3]
gene_by_cell_count_matrix <- subset(gene_by_cell_count_matrix, X..MatrixMarket <= 2500)
gene_by_cell_count_matrix <- subset(gene_by_cell_count_matrix, matrix <= 2500)
head(gene_by_cell_count_matrix)

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
colnames(data_matrix) <- cell_annotation[1:2500,1]

# Cleaning the dataset ---------------------------------------------------------

## Remove genes that are not expressed in any cell -----------------------------
keep_feature <- rowSums(data_matrix > 0) > 0
data_matrix <- data_matrix[keep_feature,]
data_matrix[1:10,1:10]

## Quality control: ------------------------------------------------------------
counts_per_cell = colSums(data_matrix)
hist(counts_per_cell, breaks = 100)

unique_genes_per_cell <- c(1:2500)

for (i in 1:ncol(data_matrix)){
  column = data_matrix[,i]
  unique_genes_per_cell[i] <- length(column[column>0]) # genes expressed
}

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

# Merge the length with the data_matrix
data_matrix <- merge(
  data_matrix,
  gene_info,
  by.x=0,
  by.y="wbps_gene_id"
)

# RPKM = numberOfReads / ( geneLength/1000 * totalNumReads/1,000,000 )
totalNumReads <- apply(data_matrix[,2:(ncol(data_matrix)-3)],2,sum)
rpkm <- matrix(, nrow = nrow(data_matrix), ncol = ncol(data_matrix)-3)
for (i in 1:nrow(data_matrix)){
  for (j in 2:(ncol(data_matrix)-3)){
    rpkm[i,j] <- 10^9 * data_matrix[i,j] / data_matrix[i,ncol(data_matrix)] / totalNumReads[j]
  }
}
rpkm <- rpkm[,2:ncol(rpkm)]
row.names(rpkm) <- data_matrix[,1]
colnames(rpkm) <- colnames(data_matrix[,2:(ncol(data_matrix)-3)])
# Generate a RPKM matrix with GO accession number and GO name at the last two columns
rpkm <- cbind(rpkm, data_matrix[,ncol(data_matrix)-2])
rpkm <- cbind(rpkm, data_matrix[,ncol(data_matrix)-1])
rpkm[1:10,1:10]

#BPSC

#Separate the cells into 5 groups based on their embryo time
cell_number<-c(1:2500)
cell_annotation<-cbind(cell_annotation,cell_number)
sorted_cell_annotation<-cell_annotation[order(cell_annotation[,"embryo.time"]),]
group1<-sorted_cell_annotation[,"cell_number"][1:500]
group2<-sorted_cell_annotation[,"cell_number"][501:1000]
group3<-sorted_cell_annotation[,"cell_number"][1001:1500]
group4<-sorted_cell_annotation[,"cell_number"][1501:2000]
group5<-sorted_cell_annotation[,"cell_number"][2001:2500]
#differential expression analysis(group1 as control group)
cell_group<-c(1:2500)
for(i in 1:ncol(rpkm)){
  if(i %in% group1){
    cell_group[i]=1
  }else if(i %in% group2){
    cell_group[i]=2
  }else if(i %in% group3){
    cell_group[i]=3
  }else if(i %in% group4){
    cell_group[i]=4
  }else{
    cell_group[i]=5
  }
}
design=model.matrix(~cell_group) 
rpkm[is.na(rpkm)] <- 0
res=BPglm(data=rpkm, controlIds=group1, design=design, coef=2, estIntPar=FALSE, useParallel=FALSE)
hist(res$PVAL, breaks=20)
#check top 5 DE genes
summary(res)
