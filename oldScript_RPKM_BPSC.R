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
