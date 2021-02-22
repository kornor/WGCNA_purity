

########### single block 

################
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 10, to=18, by=2))
# Call the network topology analysis function

sft = pickSoftThreshold(datExprPure, powerVector = powers, verbose = 3, 
                        blockSize = 5000, networkType = "signed",moreNetworkConcepts = TRUE)


# Plot the results:
sizeGrWindow(12,9)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",
     type="n",
     main = paste("Scale independence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex=cex1,
     col="red");

# this line corresponds to using an R^2 cut-off of 0.9
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",
     ylab="Mean Connectivity", 
     type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##################################################

##To check the efficiency of the chosen soft power, put it in here and plot
### Might try values as indicated by above (12 & then 14)

k=softConnectivity(datE=datExprPure,power=18, type = "signed")

# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")



## 16 = 0.9, smooth histogram (winner winner chicken dinner)
### After checking - use value of x for softthresholding
####  Proceed to module creation using automatic, dynamic cutoffs. 


softPower = 20;
adjacency = adjacency(datExprPure, power = softPower, type = "signed");
TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

######## Make the tree
geneTree = flashClust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

###### Now we're going to go ahead with module creation
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 80;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
# Calculate eigengenes
MEList = moduleEigengenes(datExprPure, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function _ no need in this case
merge = mergeCloseModules(datExprPure, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
### No merging required - all very distinct groupings
# ###################################################Rename to moduleColors (or keep dynamic colrs)
moduleColors = mergedColors
moduleColors = dynamicColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
#MEs = mergedMEs;


##### Now we'll make a dataframe to store the moduleColours for genes, 
### We'll keep adding to this baby later
### Match the names in the blocks to the colours
Modules <- as.data.frame(moduleColors)
rownames(Modules) <- colnames(datExprPure)
write.table(Modules, "Genes_moduleCall_notmerged.txt", sep = "\t")

count(Modules)
#### 

# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, Modules, 
     file = "Dynamic_pure_notmerged.RData")



##############################################################################
# load in a few things if you want to restart from here

#load(file = "Data inputs whole gene processed final.RData")
#load(file = "Dynamic_meth_merged.RData")

### kME - module membership values (kwithin)

geneModuleMembership1 = signedKME(datExprPure, MEs)
colnames(geneModuleMembership1)=paste(colnames(MEs),".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(datExprPure)[[2]]); 
colnames(MMPvalue1)=paste(colnames(MEs),".pval",sep="");


### binding the correlation and pValue tables together (optional??)
Gene       = colnames(datExprPure)
kMEtable1  = cbind(Gene,Gene,moduleColors)
for (i in 1:length(moduleColors))
  kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",
                      sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"kMEtable_notmerged_pure.csv",row.names=FALSE)
#####
write.csv(geneModuleMembership1,"MMtable__notmerged_pure.csv",row.names=TRUE)
write.csv(MMPvalue1,"Pvalues_for_MM_notmerged_pure.csv",row.names=TRUE)


##### if wish to subset by colour



black <- subset(kMEtable1, moduleColors == "black")
blue <- subset(kMEtable1, moduleColors == "blue")
brown <- subset(kMEtable1, moduleColors == "brown")
green <- subset(kMEtable1, moduleColors == "green")
greenyellow <- subset(kMEtable1, moduleColors == "greenyellow")
grey <- subset(kMEtable1, moduleColors == "grey")
purple <- subset(kMEtable1, moduleColors == "purple")
red <- subset(kMEtable1, moduleColors == "red")
turquoise <- subset(kMEtable1, moduleColors == "turquoise")

write.csv(turquoise,"kWithin_turquoise_.csv",row.names=TRUE)


#####

GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)




###########  Calculate connectivity within modules

ADJ1=abs(cor(datExprPure,use="p"))^8
Alldegrees1=intramodularConnectivity(ADJ1, moduleColors)

Alldegrees1$Modules <- moduleColors
write.table(Alldegrees1, "All_modules_connectivity_merged.txt", sep = "\t")

#######

###What about matching these modules to methylation traits???
write.table(MEs, "ModuleEigengene_Purity.txt", sep = "\t")


