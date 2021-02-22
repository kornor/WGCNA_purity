##### WGCNA Purity
## WGCNA on only high purity samples (TCGA)

###

setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_purity/WGCNA_purity")

library(WGCNA)
library(BiocManager)
library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(Hmisc)
library(blockmodeling)
library(plyr)
############################################################################################

## load the cleaned gene expression file from the previous meta-analysis
datExpr <- t(datExpr1)

## readin the purity file
pure <- read.table("High_purity_list.txt", sep = "\t", header = TRUE, row.names = 1)
## read in the normals file
normals <- read.table("Normals_list.txt", sep = "\t", header = TRUE, row.names = 1)

### trim
#### only common samples can be included
commonPure <- intersect(rownames(datExpr), rownames(pure))
datExprPure <- datExpr[commonPure,]
datExprPure <- as.data.frame(datExprPure)
commonNormal <- intersect(rownames(datExpr), rownames(normals))
datExprNormal <- datExpr[commonNormal,]
datExprNormal <- as.data.frame(datExprNormal)
###
save(datExprPure, datExprNormal, file = "Pure_Normal_trimmed_input.RData")
load(file = "Pure_Normal_trimmed_input.RData")

##########################################################################

## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(datExprPure), method = "average");


# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

## No need to remove outliers at this stage, I'm happy 
## Plot a line to show the cut




