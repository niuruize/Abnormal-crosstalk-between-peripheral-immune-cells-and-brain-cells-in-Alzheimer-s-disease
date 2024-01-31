################################################################################
###---WGCNA
################################################################################
library(impute) #BiocManager::install("impute")
library(WGCNA) #BiocManager::install("WGCNA")
library(dplyr)
library(stringr)
exp_data <- data.frame(read_xlsx("exp_data_1.xlsx"))
exp_data$median=apply(exp_data[,2:88],1,median) 
exp_data = exp_data
exp_data=exp_data[order(exp_data$Symbol,exp_data$median,decreasing = T),]
exp_data=exp_data[!duplicated(exp_data$Symbol),]
rownames(exp_data)=exp_data$Symbol
exp_data <- data.frame(exp_data[,2:88])
for (i in 1:ncol(exp_data)){exp_data[,i][is.na(exp_data[,i])] <- 0} 
expr_read <- exp_data
data.mat <- t(expr_read[order(apply(expr_read, 1, mad), decreasing = T)[1:5000],])
gsg <- goodSamplesGenes(data.mat, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(data.mat)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(data.mat)[!gsg$goodSamples], collapse = ",")));
  data.mat = data.mat[gsg$goodSamples, gsg$goodGenes]
}
#
sampleTree <- hclust(dist(data.mat), method = "average")
pdf("sampleTree_1.pdf", width = 20,height = 5)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()
#
type <- "unsigned"
powers <- c(1:10, seq(from = 12, to=30, by=2))
library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)
sft <- pickSoftThreshold(data.mat, powerVector=powers, networkType=type, verbose=3)

##
pdf("SoftPowers_1.pdf", width = 9,height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",type = "n",main = paste("Scale independence"))
text(sft$fitIndices[, 1],-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],labels = powers,cex = cex1,col = "red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[, 1],sft$fitIndices[, 5],xlab = "Soft Threshold (power)",ylab = "Mean Connectivity",type = "n",main = paste("Mean connectivity"))
text(sft$fitIndices[, 1],sft$fitIndices[, 5], labels = powers,cex = cex1,col = "red")
dev.off()
sft$powerEstimate
net <- blockwiseModules(data.mat,power = 6,maxBlockSize = 2000,TOMType = type, minModuleSize = 30,reassignThreshold = 0,
                        mergeCutHeight = 0.25,numericLabels = TRUE,pamRespectsDendro = FALSE,saveTOMs = TRUE, 
                        saveTOMFileBase = "bulk_TOM",verbose = 3)
table(net$colors)
labels2colors(net$colors)
moduleColors <- labels2colors(net$colors)
pdf("Gene dendrogram and module colors_1.pdf", width = 9,height = 5)
plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]],"Module colors", 
                    dendroLabels = FALSE,hang = 0.03,addGuide = TRUE,guideHang = 0.05)
dev.off()
MEs_col <- net$MEs
colnames(MEs_col) <- paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs_col),"ME",""))))
MEs_col <- orderMEs(MEs_col)
pdf("Eigengene adjacency heatmap_1.pdf", width = 6,height = 8)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",marDendro = c(3, 3, 2, 4),marHeatmap = c(3, 4, 2, 2),plotDendrograms = T,xLabelsAngle = 90)
dev.off()
load(net$TOMFiles[1], verbose = TRUE)
##
TOM <- as.matrix(TOM)
dissTOM <- 1-TOM
plotTOM <- dissTOM^7
diag(plotTOM) <- NA
pdf("TOMplot_1.pdf", width = 10,height = 10)
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],main = "Network heatmap plot, all genes")
dev.off()
##
mycolor <- gplots::colorpanel(250, 'red', "orange", 'lemonchiffon')
pdf("TOMplot_2.pdf", width = 10,height = 10)
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],main = "Network heatmap plot, all genes",col = mycolor)
dev.off()
##
file <- "Bulk.net"
genes <- names(net$colors[net$blockGenes[[1]]])
dimnames(TOM) <- list(genes, genes)
cyt <- exportNetworkToCytoscape(
  TOM,
  edgeFile = paste0(file, ".edges.txt"),
  nodeFile = paste0(file, ".nodes.txt"),
  weighted = TRUE,
  threshold = 0,
  nodeNames = genes,
  nodeAttr = moduleColors[net$blockGenes[[1]]])
save(net, file = "bulk.net.rda")
##
module <- "green"
moduleGenes <- names(net$colors)[which(moduleColors == module)]

##
group1 = data.frame(read_xlsx("Group.xlsx"))
use.clin = data.frame(group=group1[,2])
rownames(use.clin) = group1$sample

##
design <- model.matrix(~0 + use.clin$group)
dimnames(design) <- list(use.clin$sample, sort(unique(use.clin$group)))
design <- design[rownames(MEs_col),]

##
modTraitCor <- cor(MEs_col, design, use = "p")
modTraitP <- corPvalueStudent(modTraitCor, dim(use.clin)[1])

textMatrix <- paste0(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")")
dim(textMatrix) <- dim(modTraitCor)
pdf("clinical cor_1.pdf", width = 3,height = 7)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(design), yLabels = colnames(MEs_col),cex.lab = 0.5,ySymbols = colnames(MEs_col),
  colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)
dev.off()
##
module <- c("tan","green")
moduleGenes <- data.frame(hub_gene=names(net$colors)[which(moduleColors == module)])
write.csv(moduleGenes,"moduleGenes_hub.csv")
##
module <- c("blue")
moduleGenes <- data.frame(hub_gene=names(net$colors)[which(moduleColors == module)])
write.csv(moduleGenes,"module_blue_hub.csv")
stopCluster(cl)











