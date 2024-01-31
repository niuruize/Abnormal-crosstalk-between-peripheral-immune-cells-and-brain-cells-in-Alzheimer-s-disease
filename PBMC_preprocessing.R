################################################################################
# PBMC preprocessing
################################################################################
#rm(list=ls())
library(Seurat)
library(tidyverse)
library(liger)
library(patchwork)
library(SeuratWrappers)
library(cowplot)
library(patchwork)
library(ggplot2)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(MySeuratWrappers)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(Nebulosa) 
library(SeuratDisk)

##
data_dir <- "/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/16_对照组/filtered_feature_bc_matrix"
Control16 <- Read10X(data.dir = data_dir)
Control16 <- CreateSeuratObject(counts = Control16, project = "Control16", min.cells = 3, min.features = 200) 
Control16 <- RenameCells(Control16, add.cell.id = "Control16")
Control16[["SampleID"]] <- "Control16"
Control16[["Diagnosis"]] <- "Control"
Control16[["Age"]] <- "67"
Control16[["Sex"]] <- "F"
Control16[["percent.mt"]] <- PercentageFeatureSet(Control16,pattern = "^MT-")
summary(Control16[[]]$percent.mt)
Control16[["percent.rb"]] <- PercentageFeatureSet(Control16,pattern = "^RP[SL]")
summary(Control16[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(Control16))
Control16[["percent.HB"]] <- PercentageFeatureSet(Control16,features = HB.genes)
summary(Control16[[]]$percent.HB)
VlnPlot(Control16, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
Control16 <- subset(Control16, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 5)
VlnPlot(Control16, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#ggsave("VlnPlot.pdf", width = 28, height = 25, units = "cm")
plotControl16_1 <- FeatureScatter(Control16, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotControl16_2 <- FeatureScatter(Control16, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotControl16_1 + plotControl16_2
#ggsave("Control16_1.pdf", width = 28, height = 25, units = "cm")
rm('plotControl16_1','plotControl16_2')
Control16 <- SCTransform(Control16, vars.to.regress = "percent.mt", verbose = FALSE)
Control16 <- RunPCA(Control16)
Control16 <- FindNeighbors(Control16, dims = 1:30)
Control16 <- FindClusters(Control16, resolution = 0.8)
Control16 <- RunUMAP(Control16, dims = 1:30)
sweep.res.list <- paramSweep_v3(Control16, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Control16)*8*1e-6
homotypic.prop <- modelHomotypic(Control16$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Control16))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Control16 <- doubletFinder_v3(Control16, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Control16@meta.data)[ncol(Control16@meta.data)]="DoubletFinder"
DimPlot(Control16,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
Control16 <- CreateSeuratObject(Control16@assays$RNA@counts, meta.data = Control16@meta.data)
save(Control16,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/Control16.RData")

##
data_dir <- "/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/17_对照组/filtered_feature_bc_matrix"
Control17 <- Read10X(data.dir = data_dir)
Control17 <- CreateSeuratObject(counts = Control17, project = "Control17", min.cells = 3, min.features = 200) 
Control17 <- RenameCells(Control17, add.cell.id = "Control17")
Control17[["SampleID"]] <- "Control17"
Control17[["Diagnosis"]] <- "Control"
Control17[["Age"]] <- "73"
Control17[["Sex"]] <- "F"
Control17[["percent.mt"]] <- PercentageFeatureSet(Control17,pattern = "^MT-")
summary(Control17[[]]$percent.mt)
Control17[["percent.rb"]] <- PercentageFeatureSet(Control17,pattern = "^RP[SL]")
summary(Control17[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(Control17))
Control17[["percent.HB"]] <- PercentageFeatureSet(Control17,features = HB.genes)
summary(Control17[[]]$percent.HB)
VlnPlot(Control17, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
Control17 <- subset(Control17, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 7.5)
VlnPlot(Control17, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#ggsave("VlnPlot.pdf", width = 28, height = 25, units = "cm")
plotControl17_1 <- FeatureScatter(Control17, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotControl17_2 <- FeatureScatter(Control17, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotControl17_1 + plotControl17_2
#ggsave("Control17_1.pdf", width = 28, height = 25, units = "cm")
rm('plotControl17_1','plotControl17_2')
Control17 <- SCTransform(Control17, vars.to.regress = "percent.mt", verbose = FALSE)
Control17 <- RunPCA(Control17)
Control17 <- FindNeighbors(Control17, dims = 1:30)
Control17 <- FindClusters(Control17, resolution = 0.8)
Control17 <- RunUMAP(Control17, dims = 1:30)
sweep.res.list <- paramSweep_v3(Control17, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Control17)*8*1e-6
homotypic.prop <- modelHomotypic(Control17$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Control17))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Control17 <- doubletFinder_v3(Control17, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Control17@meta.data)[ncol(Control17@meta.data)]="DoubletFinder"
DimPlot(Control17,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
Control17 <- CreateSeuratObject(Control17@assays$RNA@counts, meta.data = Control17@meta.data)
save(Control17,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/Control17.RData")

##
data_dir <- "/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/23_早AD/filtered_feature_bc_matrix"
Early23 <- Read10X(data.dir = data_dir)
Early23 <- CreateSeuratObject(counts = Early23, project = "Early23", min.cells = 3, min.features = 200) 
Early23 <- RenameCells(Early23, add.cell.id = "Early23")
Early23[["SampleID"]] <- "Early23"
Early23[["Diagnosis"]] <- "Early_AD"
Early23[["Age"]] <- "71"
Early23[["Sex"]] <- "F"
Early23[["percent.mt"]] <- PercentageFeatureSet(Early23,pattern = "^MT-")
summary(Early23[[]]$percent.mt)
Early23[["percent.rb"]] <- PercentageFeatureSet(Early23,pattern = "^RP[SL]")
summary(Early23[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(Early23))
Early23[["percent.HB"]] <- PercentageFeatureSet(Early23,features = HB.genes)
summary(Early23[[]]$percent.HB)
VlnPlot(Early23, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
Early23 <- subset(Early23, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 7.5)
VlnPlot(Early23, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#ggsave("VlnPlot.pdf", width = 28, height = 25, units = "cm")
plotEarly23_1 <- FeatureScatter(Early23, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotEarly23_2 <- FeatureScatter(Early23, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotEarly23_1 + plotEarly23_2
#ggsave("Early23_1.pdf", width = 28, height = 25, units = "cm")
rm('plotEarly23_1','plotEarly23_2')
Early23 <- SCTransform(Early23, vars.to.regress = "percent.mt", verbose = FALSE)
Early23 <- RunPCA(Early23)
Early23 <- FindNeighbors(Early23, dims = 1:30)
Early23 <- FindClusters(Early23, resolution = 0.8)
Early23 <- RunUMAP(Early23, dims = 1:30)
sweep.res.list <- paramSweep_v3(Early23, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Early23)*8*1e-6
homotypic.prop <- modelHomotypic(Early23$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Early23))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Early23 <- doubletFinder_v3(Early23, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Early23@meta.data)[ncol(Early23@meta.data)]="DoubletFinder"
DimPlot(Early23,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
Early23 <- CreateSeuratObject(Early23@assays$RNA@counts, meta.data = Early23@meta.data)
save(Early23,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/Early23.RData")

##
data_dir <- "/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/28-早AD/filtered_feature_bc_matrix"
Early28 <- Read10X(data.dir = data_dir)
Early28 <- CreateSeuratObject(counts = Early28, project = "Early28", min.cells = 3, min.features = 200) 
Early28 <- RenameCells(Early28, add.cell.id = "Early28")
Early28[["SampleID"]] <- "Early28"
Early28[["Diagnosis"]] <- "Early_AD"
Early28[["Age"]] <- "73"
Early28[["Sex"]] <- "F"
Early28[["percent.mt"]] <- PercentageFeatureSet(Early28,pattern = "^MT-")
summary(Early28[[]]$percent.mt)
Early28[["percent.rb"]] <- PercentageFeatureSet(Early28,pattern = "^RP[SL]")
summary(Early28[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(Early28))
Early28[["percent.HB"]] <- PercentageFeatureSet(Early28,features = HB.genes)
summary(Early28[[]]$percent.HB)
VlnPlot(Early28, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
Early28 <- subset(Early28, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
VlnPlot(Early28, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#ggsave("VlnPlot.pdf", width = 28, height = 25, units = "cm")
plotEarly28_1 <- FeatureScatter(Early28, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotEarly28_2 <- FeatureScatter(Early28, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotEarly28_1 + plotEarly28_2
#ggsave("Early28_1.pdf", width = 28, height = 25, units = "cm")
rm('plotEarly28_1','plotEarly28_2')
Early28 <- SCTransform(Early28, vars.to.regress = "percent.mt", verbose = FALSE)
Early28 <- RunPCA(Early28)
Early28 <- FindNeighbors(Early28, dims = 1:30)
Early28 <- FindClusters(Early28, resolution = 0.8)
Early28 <- RunUMAP(Early28, dims = 1:30)
sweep.res.list <- paramSweep_v3(Early28, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Early28)*8*1e-6
homotypic.prop <- modelHomotypic(Early28$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Early28))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Early28 <- doubletFinder_v3(Early28, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Early28@meta.data)[ncol(Early28@meta.data)]="DoubletFinder"
DimPlot(Early28,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
Early28 <- CreateSeuratObject(Early28@assays$RNA@counts, meta.data = Early28@meta.data)
save(Early28,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/Early28.RData")

##
data_dir <- "/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/24-晚AD/filtered_feature_bc_matrix"
Late24 <- Read10X(data.dir = data_dir)
Late24 <- CreateSeuratObject(counts = Late24, project = "Late24", min.cells = 3, min.features = 200) 
Late24 <- RenameCells(Late24, add.cell.id = "Late24")
Late24[["SampleID"]] <- "Late24"
Late24[["Diagnosis"]] <- "Late_AD"
Late24[["Age"]] <- "82"
Late24[["Sex"]] <- "F"
Late24[["percent.mt"]] <- PercentageFeatureSet(Late24,pattern = "^MT-")
summary(Late24[[]]$percent.mt)
Late24[["percent.rb"]] <- PercentageFeatureSet(Late24,pattern = "^RP[SL]")
summary(Late24[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(Late24))
Late24[["percent.HB"]] <- PercentageFeatureSet(Late24,features = HB.genes)
summary(Late24[[]]$percent.HB)
VlnPlot(Late24, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
Late24 <- subset(Late24, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 7.5)
VlnPlot(Late24, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#ggsave("VlnPlot.pdf", width = 28, height = 25, units = "cm")
plotLate24_1 <- FeatureScatter(Late24, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotLate24_2 <- FeatureScatter(Late24, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotLate24_1 + plotLate24_2
#ggsave("Late24_1.pdf", width = 28, height = 25, units = "cm")
rm('plotLate24_1','plotLate24_2')
Late24 <- SCTransform(Late24, vars.to.regress = "percent.mt", verbose = FALSE)
Late24 <- RunPCA(Late24)
Late24 <- FindNeighbors(Late24, dims = 1:30)
Late24 <- FindClusters(Late24, resolution = 0.8)
Late24 <- RunUMAP(Late24, dims = 1:30)
sweep.res.list <- paramSweep_v3(Late24, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Late24)*8*1e-6
homotypic.prop <- modelHomotypic(Late24$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Late24))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Late24 <- doubletFinder_v3(Late24, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Late24@meta.data)[ncol(Late24@meta.data)]="DoubletFinder"
DimPlot(Late24,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
Late24 <- CreateSeuratObject(Late24@assays$RNA@counts, meta.data = Late24@meta.data)
save(Late24,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/Late24.RData")

##
data_dir <- "/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/27-晚AD/filtered_feature_bc_matrix"
Late27 <- Read10X(data.dir = data_dir)
Late27 <- CreateSeuratObject(counts = Late27, project = "Late27", min.cells = 3, min.features = 200) 
Late27 <- RenameCells(Late27, add.cell.id = "Late27")
Late27[["SampleID"]] <- "Late27"
Late27[["Diagnosis"]] <- "Late_AD"
Late27[["Age"]] <- "90"
Late27[["Sex"]] <- "M"
Late27[["percent.mt"]] <- PercentageFeatureSet(Late27,pattern = "^MT-")
summary(Late27[[]]$percent.mt)
Late27[["percent.rb"]] <- PercentageFeatureSet(Late27,pattern = "^RP[SL]")
summary(Late27[[]]$percent.rb)
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(Late27))
Late27[["percent.HB"]] <- PercentageFeatureSet(Late27,features = HB.genes)
summary(Late27[[]]$percent.HB)
VlnPlot(Late27, pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
Late27 <- subset(Late27, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7.5)
VlnPlot(Late27, group.by = "SampleID",pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
#ggsave("VlnPlot.pdf", width = 28, height = 25, units = "cm")
plotLate27_1 <- FeatureScatter(Late27, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotLate27_2 <- FeatureScatter(Late27, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotLate27_1 + plotLate27_2
#ggsave("Late27_1.pdf", width = 28, height = 25, units = "cm")
rm('plotLate27_1','plotLate27_2')
Late27 <- SCTransform(Late27, vars.to.regress = "percent.mt", verbose = FALSE)
Late27 <- RunPCA(Late27)
Late27 <- FindNeighbors(Late27, dims = 1:30)
Late27 <- FindClusters(Late27, resolution = 0.8)
Late27 <- RunUMAP(Late27, dims = 1:30)
sweep.res.list <- paramSweep_v3(Late27, PCs = 1:10, sct = T)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_v <- as.numeric(as.character(bcmvn$pK))
pk_good <- pk_v[bcmvn$BCmetric==max(bcmvn$BCmetric)]
DoubletRate = ncol(Late27)*8*1e-6
homotypic.prop <- modelHomotypic(Late27$seurat_clusters)
nExp_poi <- round(DoubletRate*length(colnames(Late27))) ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.agj <- round(nExp_poi*(1-homotypic.prop))
Late27 <- doubletFinder_v3(Late27, PCs = 1:10, pN = 0.25, pK = pk_good, nExp = nExp_poi.agj, reuse.pANN = FALSE, sct = T)
colnames(Late27@meta.data)[ncol(Late27@meta.data)]="DoubletFinder"
DimPlot(Late27,reduction = "umap",pt.size = 1,group.by = "DoubletFinder")
Late27 <- CreateSeuratObject(Late27@assays$RNA@counts, meta.data = Late27@meta.data)
save(Late27,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/Late27.RData")

# merge data
scRNA <- merge(Control16, y=c(Control17, Early23, Early28, Late24, Late27))
scRNA <- subset(scRNA, subset = nCount_RNA < 20000)
VlnPlot(scRNA, group.by = "SampleID", pt.size = 0,features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.HB"), ncol = 3)
plotscRNA_1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "SampleID")
plotscRNA_2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "SampleID")
plotscRNA_1 + plotscRNA_2
save(scRNA,file="/Users/niuruize/Downloads/scRNAseq/AD/AD及对照/scRNA.RData")



