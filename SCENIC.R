#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
###---SCENIC
################################################################################
library(Seurat);library(RcisTarget);library(SCENIC);library(GENIE3)
DefaultAssay(scRNA)="RNA"
table(scRNA@meta.data[["celltype"]])
colnames(scRNA@meta.data)
exprMat<-GetAssayData(object = scRNA)
dim(exprMat)
exprMat[1:4,1:10]
save(exprMat,file = "int/exprMat.Rds")
cellInfo <-data.frame(Age=scRNA$Age, celltype=scRNA$celltype,Diagnosis=scRNA$Diagnosis,cell_group=scRNA$celltype_Group)
rownames(cellInfo) <- colnames(scRNA)
head(cellInfo)
celltypeColumn <- "celltype"
colnames(cellInfo)[which(colnames(cellInfo)==celltypeColumn)] <- "celltype"
cbind(table(cellInfo$celltype))
head(cellInfo)
dir.create('int')
saveRDS(cellInfo, file="int/cellInfo.Rds")
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(celltype=c("B_Memory"="#E5D2DD",
                           "B_Naive"="#53A85F",
                           "B_Plasma"="#F1BB72",
                           "CCR6_Th"="#F3B1A0",
                           "CCR7_Th"="#D6E7A3",
                           "GZMB_Tc"="#23452F",
                           "GZMK_Tc"="#BD956A",
                           "KLRB1_Tc"="#8C549C",
                           "GZMB_KLRB1_Tc"="#585658",
                           "FOXP3+CTLA4+_Treg"="#9FA3A8",
                           "FOXP3-CTLA4+_Treg"="#E0D4CA",
                           "Platelets"="#5F3D69",
                           "NK"="#57C3F3",
                           "Monocyte"="#476D87",
                           "CMP"="#E95C59",
                           "EPC"="#E59CC4",
                           "DC"="#AB3282"))

colVars$cellcype <- colVars$celltype[intersect(names(colVars$celltype), cellInfo$celltype)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$celltype, legend=names(colVars$celltype))
##
org="hgnc" # or mgi-mouse, or dmel-fly
dbDir="/Users/niuruize/Downloads/scRNA/SCENIC/cisTarget_databases_hg19" # RcisTarget databases location
myDatasetTitle="scRNA" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
#Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
#Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
exprMat<-as.matrix(exprMat)
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat), 
                           minSamples=ncol(exprMat)*.01)
##
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
save(exprMat,file = "int/exprMat.Rds")
logMat <- log2(exprMat_filtered+1)
logMat <- na.omit(logMat)
runGenie3(logMat, scenicOptions)
##
# runSCENIC
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123
#For a very quick run:coexMethod=c("top10perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
#save...
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top10perTarget")) #** Only for toy run!!
#runSCENIC_2_createRegulons(scenicOptions) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, logMat)
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat)
savedSelections <- shiny::runApp(aucellApp)
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
scenicOptions@settings$devType="png"
runSCENIC_4_aucell_binarize(scenicOptions)
save(logMat,file = "int/logMat.Rds")
####
scenicOptions<-readRDS("int/scenicOptions.Rds")
cellInfo<-readRDS("int/cellinfo.Rds")
colVars<-readRDS("int/colVars.Rds")
exprMat<-readRDS("int/exprMat.Rds")
## group
## regulon
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_bycelltype1 <- sapply(split(rownames(cellInfo), cellInfo$celltype),
                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_bycelltype1_Scaled <- t(scale(t(regulonActivity_bycelltype1), center = T, scale=T))
pdf("celltype_TF_2.pdf", width = 5,height = 6)
pheatmap::pheatmap(regulonActivity_bycelltype1_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

# filename="regulonActivity_bycelltype1.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_bycelltype1_Scaled)
colnames(topRegulators) <- c("Regulon", "celltype", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_cell.csv")

minPerc <- .0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_bycelltype1_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$celltype),
                                                function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_bycelltype1_Binarized[which(rowSums(regulonActivity_bycelltype1_Binarized>minPerc)>0),]
pdf("celltype_TF_1.pdf", width = 5,height = 6)
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

##
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_bycelltype1 <- sapply(split(rownames(cellInfo), cellInfo$Diagnosis),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_bycelltype1_Scaled <- t(scale(t(regulonActivity_bycelltype1), center = T, scale=T))
pdf("Group_TF_1.pdf", width = 3,height = 6)
pheatmap::pheatmap(regulonActivity_bycelltype1_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()
# filename="regulonActivity_bycelltype1.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_bycelltype1_Scaled)
colnames(topRegulators) <- c("Regulon", "Group", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_group.csv") 

minPerc <- .0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_bycelltype1_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$Diagnosis),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_bycelltype1_Binarized[which(rowSums(regulonActivity_bycelltype1_Binarized>minPerc)>0),]
pdf("Group_TF_2.pdf", width = 3,height = 6)
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

#
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_bycelltype1 <- sapply(split(rownames(cellInfo), cellInfo$cell_group),
                                      function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_bycelltype1_Scaled <- t(scale(t(regulonActivity_bycelltype1), center = T, scale=T))
pdf("celltype_grpup_TF_2.pdf", width = 12,height = 7)
pheatmap::pheatmap(regulonActivity_bycelltype1_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()

# filename="regulonActivity_bycelltype1.pdf", width=10, height=20)
topRegulators <- reshape2::melt(regulonActivity_bycelltype1_Scaled)
colnames(topRegulators) <- c("Regulon", "celltype", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,"topRegulators_cell_group.csv") 

minPerc <- .0
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_bycelltype1_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$cell_group),
                                                function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_bycelltype1_Binarized[which(rowSums(regulonActivity_bycelltype1_Binarized>minPerc)>0),]
pdf("celltype_grpup_TF_1.pdf", width = 12,height = 7)
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)
dev.off()
#
rm(regulonAUC,regulonActivity_bycelltype1,regulonActivity_bycelltype1_Scaled,topRegulators,binaryRegulonActivity,
   binaryActPerc_subset,cellInfo_binarizedCells,regulonActivity_bycelltype1_Binarized)
# seurat
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(scRNA, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc.rds')

## regulonAUC
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(scRNA, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'scRNAbin.rds')

## Seurat
dir.create('scenic_seurat')
#FeaturePlot_tsne
p1 = FeaturePlot(scRNAauc, features='SPI1_670g', label=F, reduction = 'tsne')
p2 = FeaturePlot(scRNAbin, features='SPI1_670g', label=F, reduction = 'tsne')
p3 = DimPlot(scRNA, reduction = 'tsne', group.by = "Diagnosis", label=F)
p4 = DimPlot(scRNA, reduction = 'tsne', group.by = "celltype", label=F)
plotc = p1|p2|p3|p4
ggsave('scenic_seurat/SPI1_670g_tsne.png', plotc, width=19.5 ,height=4)
#FeaturePlot_umap
p1 = FeaturePlot(scRNAauc, features='SPI1_670g', label=F, reduction = 'umap')
p2 = FeaturePlot(scRNAbin, features='SPI1_670g', label=F, reduction = 'umap')
p3 = DimPlot(scRNA, reduction = 'umap', group.by = "Diagnosis", label=F)
p4 = DimPlot(scRNA, reduction = 'umap', group.by = "celltype", label=F)
plotc = p1|p2|p3|p4
ggsave('scenic_seurat/SPI1_670g_umap.png', plotc, width=19.5 ,height=4)

#RidgePlot&VlnPlot
p1 = RidgePlot(scRNAauc, features = "SPI1_670g", group.by="celltype") + theme(legend.position='none')
p2 = VlnPlot(scRNAauc, features = "SPI1_670g", pt.size = 0, group.by="Diagnosis") + theme(legend.position='none')
plotc = p1 + p2
ggsave('scenic_seurat/Ridge-Vln_SPI1_670g.png', plotc, width=10, height=8)

# pheatmap
library(pheatmap)
cellInfo <- readRDS("int/cellInfo.Rds")
celltype = subset(cellInfo,select = c('celltype','Age_stage.1'))
AUCmatrix <- t(AUCmatrix)
BINmatrix <- t(BINmatrix)
my.regulons <- rownames(AUCmatrix)
myAUCmatrix <- AUCmatrix[rownames(AUCmatrix)%in%my.regulons,]
myBINmatrix <- BINmatrix[rownames(BINmatrix)%in%my.regulons,]
# pheatmap
pheatmap(myAUCmatrix, show_colnames=F, annotation_col=celltype,
         filename = 'scenic_seurat/myAUCmatrix_heatmap_1.pdf',
         width = 5, height = 8)
# pheatmap
pheatmap(myBINmatrix, show_colnames=F, annotation_col=celltype,
         filename = 'scenic_seurat/myBINmatrix_heatmap.pdf',
         color = colorRampPalette(colors = c("white","black"))(100),
         width = 5, height = 6)
rm(p1,p2,p3,p4,myAUCmatrix,myBINmatrix,scRNA,scRNAauc,scRNAbin,AUCmatrix,BINmatrix,RegulonName_AUC,RegulonName_BIN,
   regulonActivity_bycelltype1,regulonActivity_bycelltype1_Binarized,regulonActivity_bycelltype1)











