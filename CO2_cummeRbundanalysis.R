# This is an R script to analyze 
# the RNA-Seq data from an analysis 
# of the carbon concentrating mechanism
# in Chlamydomonas reinhardtii as
# described in Fang et al., TPC 2012
# www.plantcell.org/cgi/doi/10.1105/tpc.112.097949

# For this analysis, we will be using
# a collection of functions 
# (i.e. a library) of tools specifically
# for analyzing RNA-Seq data called cummeRbund.
# Load the cummeRbund library with this:
library(cummeRbund)

# For this script, we will also need 
# the dplyr library:
library(dplyr)

# Next we are going to set the working directory with setwd()
setwd("~/MCB_C117_Lecture_21/")

# Next, we need to load the RNA-Seq data.
cuff<-readCufflinks(dir = "RNAseq_data_CO2/")


# To check that the data has been processed 
# correctly, we can type the name of the 
# database and it should return a brief
# summary of the data:

cuff

# You should see
# CuffSet instance with:
#	 6 samples 
#	 17741 genes etc...


# Make a table of FPKMs for all genes
fpkms<-fpkmMatrix(genes(cuff))

# Check the first few lines
head(fpkms)

# Export the table of FPKMs
geneIDs<-row.names(fpkms)

write.table(cbind(geneIDs, fpkms), file = "CO2_FPKMs.tsv", sep = "\t", row.names=FALSE)


# Check expression of the CIA5 gene.
# The gene ID for CIA5 is Cre02.g096300.
# Put the gene ID into a character variable

cia5<-"Cre02.g096300"

# Select the cia5 row in our fpkms data frame
fpkms[cia5,]


# plot a dendrogram comparing samples
csDendro(genes(cuff))

# plot a dendrogram with separate replicates
csDendro(genes(cuff),replicates=TRUE)


# Make a PCA plot of each condition (called MDS in cummeRbund)
MDSplot(genes(cuff),replicates = FALSE)

# Make a PCA plot of each replicate
MDSplot(genes(cuff),replicates = TRUE)

# examine differential expression data
gene.diff<-diffData(genes(cuff))

# check the first few lines
head(gene.diff)

# how many pairwise comparisons are there?
dim(gene.diff)
# 17,741 genes x 15 comparisons = 266,115!

# get all genes that are 
# significantly different in any sample
mySigGeneIds<-getSig(cuff,level='genes')

# check the first few lines
head(mySigGeneIds)

# how many genes are significant?
length(mySigGeneIds)

# get significant genes just
# between WT_VLCO2 vs WT_HiCO2
WTsigGenes<-getSig(cuff,
                  x='WT_VLCO2',
                  y='WT_HiCO2',
                  level='genes')

# how many genes in that comparison?
length(WTsigGenes)

# get all gene data for those genes
WTsigGeneData<-getGenes(cuff,WTsigGenes)

# make a basic heatmap of those genes with default settings
csHeatmap(WTsigGeneData)

# tweak the settings for the heatmap
csHeatmap(WTsigGeneData,
             labRow = FALSE,
             heatscale = c("yellow","orange","red"),
             logMode = TRUE)

csHeatmap(WTsigGeneData,
             labRow = FALSE,
             heatscale = c("gray90","darkblue"),
             logMode = TRUE)

# Get significant genes just
# between cia5_VLCO2 vs WT_VLCO2
cia5vsWTsigGenes<-getSig( cuff,
                          x='cia5_VLCO2',
                          y='WT_VLCO2',
                          level='genes')

# How many genes in that comparison?
length(cia5vsWTsigGenes)
# How does that compare to the previous comparison?

# Perform k-means clustering with k=4 on WTsigGenes
cluster_k4<-csCluster(WTsigGeneData,k=4)

# Take a look at the results
head(cluster_k4$cluster)

# Make a plot of the clusters 
csClusterPlot(cluster_k4)

# get cluster 4 gene IDs
myClusters<-as.data.frame(cluster_k4$clustering)
cluster4<-row.names(filter(myClusters, cluster_k4$cluster == 4))


# how many genes are in cluster 4?
length(cluster4)

# export the cluster4 gene IDs
write(cluster4, file = "cluster_4_gene_IDs.txt")






