library("BPSC")
library("SingleCellExperiment")
library("scater")
library("biomaRt")
library("Matrix")
library('scuttle')
library("pheatmap")
library("MAST")
library("iDEA")

load("workspace.RData")

gene_by_cell_count_matrix <- read.csv("GSE126954_gene_by_cell_count_matrix.txt", sep = "")[,1:3]

gene_annotation <- read.csv("GSE126954_gene_annotation.csv")[,2:3]
cell_annotation <- read.csv("GSE126954_cell_annotation.csv")[,-1]

data_matrix = sparseMatrix(i=gene_by_cell_count_matrix$X..MatrixMarket,
                          j=gene_by_cell_count_matrix$matrix,
                          x=gene_by_cell_count_matrix$coordinate
                           )
rownames(data_matrix)=gene_annotation[,1]
colnames(data_matrix)=cell_annotation[,1]

head(data_matrix)
dim(data_matrix)

cells_bin = cell_annotation[cell_annotation$raw.embryo.time.bin=="330-390",]
data_matrix = data_matrix[,as.numeric(unlist(rownames(cells_bin)))]

sce <- SingleCellExperiment(
  assays = list(counts = data_matrix),
  rowData = data.frame(gene_names = rownames(data_matrix)),
  colData = data.frame(cell_names = colnames(data_matrix))
)

## Remove genes that are not expressed in any cell
keep_feature <- rowSums(counts(sce) > 0) > 0
sce_filtered <- SingleCellExperiment(
  assays = list(counts = counts(sce)[keep_feature, ]),
  rowData = data.frame(gene_names = rowData(sce)[keep_feature,]),
  colData = data.frame(cell_names = colData(sce))
)

sce_filtered <- addPerCellQC(sce_filtered)
sce_filtered <- addPerFeatureQC(sce_filtered)
hist(sce_filtered$total, breaks = 100) #hist of counts per cell

unique_genes_per_cell <- colSums(counts(sce)>0)
hist(unique_genes_per_cell,breaks=100) #hist of unique genes per cell

cpm(sce_filtered) <- calculateCPM(sce_filtered) # using scater package
sce_filtered <- logNormCounts(sce_filtered,transform="log") 
# this logNormCounts function comes from scuttle package. Do not change sce_filtered into logcounts(sce_filtered) although that seems more logical.
# The CPM and log counts are now stored in this sce_filtered SCE object
# str(sce_filtered) shows the hierarchical structure of this SCE object

# Sparse matrices for further analysis
sce_cpm <- cpm(sce_filtered) # CPM for BPSC input
sce_log <- sce_filtered@assays@data@listData$logcounts # log transformed counts for heatmap input

cell_type <- cell_annotation[cell_annotation$cell %in% dimnames(sce_cpm)[[2]],6] 

unique_cell_types <- unique(cell_type)[!is.na(unique(cell_type))]
unique_cell_types
length(unique_cell_types)

cell_indices <- list()
# Extract the cell indices according to each of the cell type and put them in the list
for (i in 1:length(unique_cell_types)){
  cell_indices[[unique_cell_types[i]]] <- dimnames(sce_filtered[,which(cell_type==unique_cell_types[i])])[[2]]
}

cell_indices[[1]][1:5]
cell_indices[['Body_wall_muscle']][1:5]

matList <- list()
cellTypeList <- list()

for (i in 1:length(unique_cell_types)){
  mat <- sce_log[1:20,(cell_indices[[i]])] #input:logcounts
  matList[[i]] <- mat
  if(is.null(ncol(mat))){
    cellTypeList[[i]] <- unique_cell_types[i]
  } else {
    cellTypeList[[i]]<-rep(unique_cell_types[i],ncol(mat)) 
  }
}

mat.whole <- do.call(cbind,matList)
cellType <- do.call(c,cellTypeList)
cellType <- data.frame(cellType)
dim(mat.whole)[2] == nrow(cellType)
colnames(cellType) <- "Cell type"
mat.whole <- data.frame(as.matrix(mat.whole))
rownames(cellType) <- colnames(mat.whole)

pheatmap(mat.whole,annotation=cellType,cluster_cols=FALSE)

control.countData=sce_filtered[1:20,cell_indices[[1]]]
treated.countData=sce_filtered[1:20,cell_indices[[2]]]
countData <- cbind(control.countData,treated.countData)
cellType <- rep(c(1,2), c(ncol(control.countData),ncol(treated.countData)))
cellType <- as.data.frame(cellType)
colnames(cellType) <- "CellType"
colData(countData) <- cbind(colData(countData),cellType)
head(colData(countData),3)

colData(countData)$CellType<-factor(colData(countData)$CellType)
countData<-SceToSingleCellAssay(countData, class = "SingleCellAssay",check_sanity = FALSE) #uses MAST package
zlmCond <- zlm(~CellType, countData)
logFC<-getLogFC(zlmCond)[,c(1,3)]

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

# Use biomRT package to get access to WormBase ParaSite BioMart
# Establish a connection to the WormBase ParaSite BioMart
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# Retrieve the GO terms
gene_info <- getBM(filters="wbps_gene_id", 
  attributes=c("wbps_gene_id", "go_accession","go_name_1006"), 
  values=row.names(data_matrix), 
  mart=mart,
  uniqueRows=FALSE)

dim(gene_info) 

head(gene_info)

#Here the first results(only two cell types) is used!!
pvalue <- results$`P-value` #### the pvalue column
zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
fc <- results$logFC ## the fold change column
se_beta <- abs(fc/zscore) ## to approximate the standard error of beta
var = se_beta^2  ### square 

# Summary is a matrix of fold change and variance of each gene
summary = data.frame(fc = fc,variance = var,row.names = results$gene)

# get annotation data
# know how many go terms we have
length(unique(gene_info$go_name))#6907
if ("annotation.RData" %in% dir()){
    load(file="annotation.RData")
} else {
    annotation<-matrix(0,nrow =nrow(summary) ,ncol = 6907)
    rownames(annotation)<-results$gene
    colnames(annotation)<-unique(gene_info$go_name_1006)
    for (i in 1:nrow(annotation)) {
      for (j in 1:ncol(annotation)) {
        index<-which(gene_info$wbps_gene_id==rownames(annotation)[i])
        for (k in 1:length(index)) {
          if(gene_info$go_name_1006[k]==colnames(annotation)[j]) {
            annotation[i,j]=1
          }
        }
      }
    }
    save(annotation, file="annotation.RData")
}

#create idea object
idea<-CreateiDEAObject(summary,annotation)

#Fit the model
options(mc.cores = 1) # to fix an error. Has to do with problems in parallelization 
idea <- iDEA.fit(idea)

#correct p-values
idea <- iDEA.louis(idea)

#get output
idea@gsea

#DE analysis of individual gene
#with pre-selected genes
pip = unlist(idea@de[["membrane"]]$pip)
head(pip)

#without pre-selected genes
idea <- iDEA.BMA(idea)
head(idea@BMA_pip)

save.image(file="workspace.RData")

genes_of_interest = as.matrix(read.delim("genes_of_interest.txt", head=FALSE))
genes_index = which(gene_annotation[,2] %in% genes_of_interest)
geneIDs = gene_annotation[genes_index,1]

combinations = data.frame()
for (i in 1:(length(unique_cell_types)-1)){
    for (j in (i+1):length(unique_cell_types)){
        row = c(i,j)
        combinations = rbind(combinations, row)
    }
}

annotation = sparseMatrix(i=nrow(sce_cpm),
                          j=length(unique(gene_info$go_name)),
                          x=0
                           )
rownames(annotation)<-rownames(sce_cpm)
colnames(annotation)<-unique(gene_info$go_name_1006)

for (i in 1:nrow(annotation)) {
    if (i%%100==0) i
    index<-which(gene_info$wbps_gene_id==rownames(annotation)[i])
    if (length(index)==0) next
    for (j in 1:ncol(annotation)) {
      for (k in 1:length(index)) {
        if(gene_info$go_name_1006[k]==colnames(annotation)[j]) {
          annotation[i,j]=1
        }
      }
    }
}

annotation_matrix = as.matrix(annotation)
annotation_df = as.data.frame(annotation_matrix)
annotation = annotation_df

# Create a list for log fold change + p-value from BPSC results and for iDEA results
results <- list()
list_of_idea = list()

options(mc.cores = 1) # to fix an error. Has to do with problems in parallelization in idea

for (i in 1:nrow(combinations)){
  print(i)
  startTime = Sys.time()
  #Define the two groups to be compared (Remember of remove 1:20 when we run it for the whole data containing all the cells in a time bin!!!)
  control.mat=sce_cpm[geneIDs,cell_indices[[combinations[i,1]]]]
  treated.mat=sce_cpm[geneIDs,cell_indices[[combinations[i,2]]]]
  #Create a data set by merging the control group and the treated group
  bp.mat=cbind(control.mat,treated.mat)
  rownames(bp.mat)=c(1:nrow(bp.mat))
  colnames(bp.mat)=c(1:ncol(bp.mat))
  group=c(rep(1,ncol(control.mat)),rep(2,ncol(treated.mat)))
  #Run BPglm for differential expression analysis
  res=BPglm(data=bp.mat, controlIds=which(lapply(group, as.numeric)==1), design=model.matrix(~group), coef=2, estIntPar=FALSE, useParallel=FALSE)
  
  # Use MAST to compute logFC based on the raw counts
  control.countData=sce_filtered[geneIDs,cell_indices[[combinations[i,1]]]]
  treated.countData=sce_filtered[geneIDs,cell_indices[[combinations[i,2]]]]
  countData <- cbind(control.countData,treated.countData)
  cellType <- as.data.frame(rep(c(1,2), c(ncol(control.countData),ncol(treated.countData))))
  colnames(cellType) <- "CellType"
  colData(countData) <- cbind(colData(countData),cellType)
  colData(countData)$CellType<-factor(colData(countData)$CellType)
  countData<-SceToSingleCellAssay(countData, class = "SingleCellAssay",check_sanity = FALSE)
  zlmCond <- zlm(~CellType, countData)
  logFC<-getLogFC(zlmCond)[,c(1,3)]
    
  # Log fold change + p-value from BPSC
  PVAL<-as.data.frame(res$PVAL)
  result <- cbind(logFC,res$PVAL)
  colnames(result) <- c('gene','logFC','P-value')
  results[[paste(unique_cell_types[combinations[i,1]],unique_cell_types[combinations[i,2]], sep=" vs. ")]]=result
  
  ###iDEA
  #calculate variance
  #Here the first results(only two cell types) is used!!
  pvalue <- result$`P-value` #### the pvalue column
  zscore <- qnorm(pvalue/2.0, lower.tail=FALSE) #### convert the pvalue to z-score
  fc <- result$logFC ## the fold change column
  se_beta <- abs(fc/zscore) ## to approximate the standard error of beta
  var = se_beta^2  ### square 
  summary = data.frame(fc = fc,variance = var,row.names = result$gene)# Summary is a matrix of fold change and variance of each gene
  
  #create idea object
  idea<-CreateiDEAObject(summary,annotation[geneIDs,])
  #Fit the model
  idea <- iDEA.fit(idea)
  #correct p-values
  idea <- iDEA.louis(idea)
  #DE analysis of individual gene
  #with pre-selected genes
  pip = unlist(idea@de[["membrane"]]$pip)
  #without pre-selected genes
  idea <- iDEA.BMA(idea)
    
  #Save idea to list, so we can do analysis later
  list_of_idea <- cbind(list_of_idea,idea)
  
  stopTime = Sys.time()
  print(stopTime-startTime)
} 
