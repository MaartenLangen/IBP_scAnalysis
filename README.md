# IBP_scAnalysis
Repository for all code related to our IBP project. The jupyter notebooks contain the complete analysis and the code for making the visualizations.

# Included files
The cell and gene annotation files are included. The file with counts for the genes in each cell is too large to put on github and should thus be downloaded independently. [This is a link to the count matrix.](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE126954&format=file&file=GSE126954%5Fgene%5Fby%5Fcell%5Fcount%5Fmatrix%2Etxt%2Egz)
A file with genes of interest is also supplied: it is necessary for part of our analysis and contains a list of genes that might be related to the cytoskeleton.

# Running the jupyter notebooks
All commands used to create an environment in which the jupyter notebooks with R can be ran, are in the creatingEnvironment.txt file. The packages that should be present in this environment, can be installed using the packegesInEnvironment.yml file.

# More on the notebooks
GSEA_bin330.ipynb contains code for the Gene Set Enrichment Analysis performed in this project. It uses the data of timebin 330-390 for this, but this can be easily changed to include other timebins. It generates a count matrix based on the data by [Packer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE126954). The notebook ends with creating a bubble plot.
Heatmap_spatial.ipynb contains the code for differential gene expression in a particular time bin for different comparisons between cells. Heatmap_temporal.ipynb is similar to the previous notebook, but now looks at particular cell types and visualizes expression over time.
