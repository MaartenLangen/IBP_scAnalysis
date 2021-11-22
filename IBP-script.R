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
cell_type <- subset(cell_annotation, cell_annotation$cell %in% colnames(sce_cpm))[,6]

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
res=BPglm(data=bp.mat, controlIds=controlIds, design=design, coef=coef, estIntPar=FALSE, useParallel=FALSE, keepFit = FALSE) 
#Plot the p-value distribution
res$PVAL
hist(res$PVAL, breaks=20)
#Summarize the resutls
summary(res)
#Fold change
fc=apply(treated.mat,1,mean)/apply(control.mat,1,mean)
fc # the NA values indicates division by zero
# BPSC + fold change results
results <- t(rbind(fc,res$PVAL))
colnames(results) <- c('Fold change','P-value')
results

# use BPSC to calculate fold change and variance
# To do this, first we need to fit parameters to all sets (independent for control and treated)
mat.control.res = estimateBPMatrix(dataMat = control.mat, para.num = 4, estIntPar = TRUE)
mat.treated.res = estimateBPMatrix(dataMat = treated.mat, para.num = 4, estIntPar = TRUE)
# put the parameters in nice matrices
pc = matrix(unlist(mat.control.res$bp.model.list[names(mat.control.res$bp.model.list)=="par"]), ncol=4, byrow=TRUE)
pt = matrix(unlist(mat.treated.res$bp.model.list[names(mat.treated.res$bp.model.list)=="par"]), ncol=4, byrow=TRUE)
k=1
l=1
param.control = NULL
param.treated = NULL
for (i in 1:nrow(control.mat)){
  if (i %in% mat.control.res$ind.set){
    param.control = rbind(param.control, pc[k,])
    k = k+1
  } else {
    param.control = rbind(param.control, c(NA, NA, NA, NA))
  }
  if (i %in% mat.treated.res$ind.set){
    param.treated = rbind(param.treated, pt[l,])
    l = l+1
  } else {
    param.treated = rbind(param.treated, c(NA, NA, NA, NA))
  }
}
rownames(param.control) = rownames(control.mat)
rownames(param.treated) = rownames(treated.mat)

# use these matrices to calculate fold change and variance
fc = NULL
fc.var = NULL
for (i in 1:nrow(param.control)){
  mean.control = meanBP(alp=param.control[i,1],bet=param.control[i,2],lam1=param.control[i,3], lam2=param.control[i,4])
  mean.treated = meanBP(alp=param.treated[i,1],bet=param.treated[i,2],lam1=param.treated[i,3], lam2=param.treated[i,4])
  var.control = varBP(alp=param.control[i,1],bet=param.control[i,2],lam1=param.control[i,3], lam2=param.control[i,4])
  var.treated = varBP(alp=param.treated[i,1],bet=param.treated[i,2],lam1=param.treated[i,3], lam2=param.treated[i,4])
  foldChange = mean.control/mean.treated
  varFoldChange = foldChange*sqrt((sqrt(var.control)/mean.control)^2+(sqrt(var.treated)/mean.treated)^2)
  fc=rbind(fc,c(foldChange,varFoldChange))
}
colnames(fc) = c("foldChange","variance")


## Generalize the small test to all possible combinations of cell types
# Generate all combinations of numbers from 1 to the number of cell types, taken 2 at a time
combinations <- expand.grid(1:length(unique_cell_types), 1:length(unique_cell_types))
head(combinations)
combinations <- combinations[combinations$Var1 != combinations$Var2,]
head(combinations)
nrow((combinations)) # number of combinations
# Create a list for BPSC + fold change results
results <- list()
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
  #Fold change
  fc=apply(treated.mat,1,mean)/apply(control.mat,1,mean)
  fc
  # BPSC + fold change results
  result <- t(rbind(fc,res$PVAL))
  colnames(result) <- c('Fold change','P-value')
  results[[paste(unique_cell_types[combinations[i,1]],unique_cell_types[combinations[i,2]], sep=" vs. ")]]=result
} 



### Code that may be useful later on
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

# Calculate the gene length (do not need it)
gene_info$length <- (gene_info$end_position - gene_info$start_position)
gene_info <- gene_info[-c(2:3)]
colnames(gene_info)[3] <- "go_name"
head(gene_info)