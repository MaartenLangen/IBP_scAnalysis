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
library("biomaRt")
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

# select a small part of data (the first 1000 cells in the matrix within the selected time bin)
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

# Normalize data (CPM) using scater package
sce_cpm <- calculateCPM(sce_filtered)
sce_cpm[1:5,1:5]

### Generate the column indices of the gene expression matrix corresponding to each cell type 

# First, generate the cell types included in the matrix being analyzed 
# and make them correspond to the cell names (column names of the matrix)
# cell_annotation[,6] contains the cell types we need.
# We use them to split the data into groups containing gene expression of different cell types
# dimnames(sce_cpm)[[2]] contains the column names of the sparse matrix
cell_type <- cell_annotation[cell_annotation$cell==dimnames(sce_cpm)[[2]],6] 

# Unique values in the vector cell_type that is not NA
# We do not want to take the type-unidentified (na) cells into the differential expression analysis
unique_cell_types <- unique(cell_type)[!is.na(unique(cell_type))] # we do not take the type-unidentified (na) cells into the differential expression analysis
unique_cell_types
length(unique_cell_types)

# Create a list for cell indices
cell_indices <- list()
# Extract the cell indices according to each of the cell type and put them in the list
for (i in 1:length(unique_cell_types)){
  cell_indices[[unique_cell_types[i]]] <- dimnames(sce_cpm[,which(cell_type==unique_cell_types[i])])[[2]]
}
# We can retrieve the cell indices of a cell type  by either the index of the list, or the cell type such as 'Body_wall_muscle'
cell_indices[[1]][1:5]
cell_indices[['Body_wall_muscle']][1:5]


### Use MAST to compute logFC based on the raw counts
control.countData=sce_filtered[1:20,cell_indices[[1]]]
treated.countData=sce_filtered[1:20,cell_indices[[2]]]
countData <- cbind(control.countData,treated.countData)
cellType <- rep(c(1,2), c(ncol(control.countData),ncol(treated.countData)))
cellType <- as.data.frame(cellType)
colnames(cellType) <- "CellType"
colData(countData) <- cbind(colData(countData),cellType)
head(colData(countData),3)
colData(countData)$CellType<-factor(colData(countData)$CellType)
countData<-SceToSingleCellAssay(countData, class = "SingleCellAssay",check_sanity = FALSE)
zlmCond <- zlm(~CellType, countData)
logFC<-getLogFC(zlmCond)[,c(1,3)]
logFC

### Differential expression analysis using BPSC 
# (Here we take only the first 20 genes to run a quick small test because BPSC takes a bit long time to run)

## A test to compare the gene expression of only the first two cell types: "Body_wall_muscle" and "Ciliated_amphid_neuron"
#Define the two groups to be compared
control.mat=sce_cpm[1:20,cell_indices[[1]]]
treated.mat=sce_cpm[1:20,cell_indices[[2]]]
#Create a data set by merging the control group and the treated group
bp.mat=cbind(control.mat,treated.mat)
rownames(bp.mat)=c(1:nrow(bp.mat))
colnames(bp.mat)=c(1:ncol(bp.mat))
group=c(rep(1,ncol(control.mat)),rep(2,ncol(treated.mat)))
#First, choose IDs of all cells of the control group for estimating parameters of BP models
controlIds=which(group==1)
#Create a design matrix including the group labels. All batch effects can be also added here if they are available
design=model.matrix(~group) 
#Select the column in the design matrix corresponding to the coefficient (the group label) for the GLM model testing
coef=2 
#Run BPglm for differential expression analysis
res=BPglm(data=bp.mat, controlIds=controlIds, design=design, coef=coef, estIntPar=FALSE, useParallel=FALSE) 
#Plot the p-value distribution
res$PVAL
hist(res$PVAL, breaks=20)
#Summarize the resutls
summary(res)

# Log fold change + p-value from BPSC results
PVAL<-as.data.frame(res$PVAL)
results <- cbind(logFC,res$PVAL)
colnames(results) <- c('gene','logFC','P-value')
results



### Code that may be useful later on
# Use biomRT package to get access to WormBase ParaSite BioMart
# Establish a connection to the WormBase ParaSite BioMart
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# List available filters that can be used in a query and attributes you can retrieve
listFilters(mart)[1:100,]
listAttributes(mart)[1:100,]

# Retrieve the GO terms
gene_info <- getBM(filters="wbps_gene_id", 
  attributes=c("wbps_gene_id", "go_accession","go_name_1006"), 
  values=row.names(data_matrix), 
  mart=mart,
  uniqueRows=FALSE)
dim(gene_info) 
head(gene_info)




###iDEA
#calculate variance
#Here the first results(only two cell types) is used!!
pvalue <- results[,2] #### the pvalue column
zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
fc <- results[,1] ## the fold change column
se_beta <- abs(fc/zscore) ## to approximate the standard error of beta
var = se_beta^2  ### square 
summary = data.frame(fc = fc,variance = var)# Summary is a matrix of fold change and variance of each gene
#get annotation data
#know how many go terms we have
length(unique(gene_info$go_name))#2768
annotation<-matrix(0,nrow =nrow(summary) ,ncol = 2768)
rownames(annotation)<-rownames(summary)
colnames(annotation)<-unique(gene_info$go_name)
for (i in 1:nrow(annotation)) {
  for (j in 1:ncol(annotation)) {
    index<-which(gene_info$wbps_gene_id==rownames(annotation)[i])
    for (k in 1:length(index)) {
      if(gene_info$go_name[k]==colnames(annotation)[j]) {
        annotation[i,j]=1
    }
    
    }
  }
 
}

#install iDEA package
devtools::install_github('xzhoulab/iDEA')
library(iDEA)
#create idea object
idea<-CreateiDEAObject(summary,annotation)
#Fit the model
idea <- iDEA.fit(idea)
#correct p-values
idea <- iDEA.louis(idea)
#get output
idea@gsea


## Generalize the small test to all possible combinations of cell types
# Generate all combinations of numbers from 1 to the number of cell types, taken 2 at a time
combinations <- expand.grid(1:length(unique_cell_types), 1:length(unique_cell_types))
head(combinations)
combinations <- combinations[combinations$Var1 != combinations$Var2,]
head(combinations)
nrow((combinations)) # number of combinations

#Until now it's just for comparison of two cell types, so the following code is for adding the iDea analysis to the "for loop" for several pair-wise comparisons.
# Create a list for log fold change + p-value from BPSC results
results <- list()
results_iDEA<-list()
for (i in 1:nrow(combinations)){
  #Define the two groups to be compared (Remember of remove 1:20 when we run it for the whole data containing all the cells in a time bin!!!)
  control.mat=sce_cpm[1:20,cell_indices[[combinations[i,1]]]]
  treated.mat=sce_cpm[1:20,cell_indices[[combinations[i,2]]]]
  #Create a data set by merging the control group and the treated group
  bp.mat=cbind(control.mat,treated.mat)
  rownames(bp.mat)=c(1:nrow(bp.mat))
  colnames(bp.mat)=c(1:ncol(bp.mat))
  group=c(rep(1,ncol(control.mat)),rep(2,ncol(treated.mat)))
  #Run BPglm for differential expression analysis
  res=BPglm(data=bp.mat, controlIds=which(lapply(group, as.numeric)==1), design=model.matrix(~group), coef=2, estIntPar=FALSE, useParallel=FALSE)
  
  # Use MAST to compute logFC based on the raw counts
  control.countData=sce_filtered[1:20,cell_indices[[combinations[i,1]]]]
  treated.countData=sce_filtered[1:20,cell_indices[[combinations[i,2]]]]
  countData <- cbind(control.countData,treated.countData)
  cellType <- as.data.frame(rep(c(1,2), c(ncol(control.countData),ncol(treated.countData))))
  colnames(cellType) <- "CellType"
  colData(countData) <- cbind(colData(countData),cellType)
  colData(countData)$CellType<-factor(colData(countData)$CellType)
  countData<-SceToSingleCellAssay(countData, class = "SingleCellAssay",check_sanity = FALSE)
  zlmCond <- zlm(~cellType, countData)
  logFC<-getLogFC(zlmCond)[,c(1,3)]
  
  # Log fold change + p-value from BPSC
  PVAL<-as.data.frame(res$PVAL)
  result <- cbind(logFC,res$PVAL)
  colnames(results) <- c('gene','logFC','P-value')
  results[[paste(unique_cell_types[combinations[i,1]],unique_cell_types[combinations[i,2]], sep=" vs. ")]]=result
  
  #iDEA
  #calculate variance
  pvalue <- res$PVA
  zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
  se_beta <- abs(fc/zscore) ## to approximate the standard error of beta
  var = se_beta^2  ### square 
  summary = data.frame(fc = fc,variance = var)# Summary is a matrix of fold change and variance of each gene
  #create idea object
  idea<-CreateiDEAObject(summary,annotation)
  #Fit the model
  idea <- iDEA.fit(idea)
  #correct p-values
  idea <- iDEA.louis(idea)
  #Save all results matrix in a list
  results_iDEA[[i]]<-idea@gsea
} 






