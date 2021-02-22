
## normals take two!

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

# load the PURITY WGCNA modules

load(file = "Dynamic_pure_notmerged.RData")
load(file = "Pure_Normal_trimmed_input.RData")

######### Calculation of module preservation
datExprPure <- as.data.frame(t(datExprPure))
datExprNorm <- as.data.frame(t(datExprNorm))


### preprocess the data using goodgenes (WCGNA)

gsg = goodSamplesGenes(datExprNorm,verbose = 5);
gsg$allOK

## if return is true, all good to go
### Otherwise

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprNorm)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprNorm)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprNorm = datExprNorm[gsg$goodSamples, gsg$goodGenes]
}

### and for the pures

gsg = goodSamplesGenes(datExprPure,verbose = 5);
gsg$allOK

if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExprPure)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExprPure)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprPure = datExprPure[gsg$goodSamples, gsg$goodGenes]
}


########## check that gene names match 
##Pure now the active expression file
commongenes <- intersect(colnames(datExprPure), colnames(datExprNorm))
datExprPure <- datExprPure[,commongenes]

##trim the module colours
Module1 <- Modules[commongenes,]


setLabels = c("Pure","Norm");
multiExpr = list(Pure= list(data = datExprPure), Norm = list(data = datExprNorm));
multiColor = list(Pure = Module1);

### check for preservation

system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation_Normals_vs_Pure.RData");
load(file = "modulePreservation_Normals_vs_Pure.RData")

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

# Compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

pqual <- as.data.frame(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)))

write.table(pqual, "Preservation vs Quality for Norm to Pure.txt", sep = "\t")

##graphing
# Module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];

# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
# Main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/TCGA_PDX-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();


#### now to calculate the module eigengene for each module as in TCGA:

## first need to give mod_genes the row names matching genes of TCGA set
row.names(mod_gene) <- colnames(datExprNorm)

# now trim to just the ones in the mod_normal set

commongenes <- intersect(rownames(mod_gene), colnames(datExprNorm))

mod_normal <- mod_gene[commongenes,]
#mod_normal <- as.data.frame(mod_normal)
#rownames(mod_normal) <- commongenes


###now the calcululation:

PCs_normal <- moduleEigengenes(datExprNorm, 
                               mod_normal, 
                               impute = TRUE, 
                               nPC = 1, 
                               align = "along average", 
                               excludeGrey = FALSE, 
                               grey = if (is.numeric(mod_pdx)) 0 else "grey",
                               subHubs = TRUE,
                               softPower = 6,
                               scale = TRUE,
                               verbose = 0, indent = 0)

####works!!

ME_normals    = PCs_normal$eigengenes
distPC_normals = 1-abs(cor(ME_normals,use="p"))
distPC_normals = ifelse(is.na(distPC_normals), 0, distPC_normals)
pcTree_normals = hclust(as.dist(distPC_normals),method="a") 
MDS_normals  = cmdscale(as.dist(distPC_normals),2)
colorsNormals = names(table(mod_normal))



### add in the rownames

#rownames(ME_normals) <- colnames(datExprNorm)

write.table(ME_normals, file = "ME_Normals.txt", sep = "\t")
save(PCs_normal, file = 'Normals_module_preservation_data.Rdata')




##### MDS of the module expression in the PDX group:
plot.new()
# Start the plot
MDS_N1 <-as.data.frame(MDS_normals)
sizeGrWindow(12, 9);
plot(MDS_N1, col= colorsNormals, main="MDS plot Normals", cex=2, pch=19, xlab = "Principal component 1", ylab = "Principal component 2")

plot(MDS_normals$V1, MDS_normals$V2)

######### look at the samples clustering?

## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(brca_t), method = "average");


# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)




##############################
greens <- read.table("Green_sets.txt", sep = "\t", header = TRUE)


greens$Pure_set_Green %in% greens$TCGA_set_Green

a <- intersect(greens$Pure_set_Green,greens$TCGA_set_Green) # in both
b <- setdiff(greens$Pure_set_Green,greens$TCGA_set_Green) # in A but not B
c <- setdiff(greens$TCGA_set_Green,greens$Pure_set_Green) # in B but not A


write.csv(a, file = "Joint_greens.csv")
write.csv(b, file = "Pure_greens.csv")
write.csv(c, file = "TCGA_greens.csv")
