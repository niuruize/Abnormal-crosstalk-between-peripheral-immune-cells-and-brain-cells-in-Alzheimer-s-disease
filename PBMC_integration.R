# loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
# PBMC data analysis
################################################################################
load("scRNA.RData")
scRNA = subset(scRNA, DoubletFinder == "Singlet")
{
  scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
  scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
  })
  scRNA.features <- SelectIntegrationFeatures(object.list = scRNAlist, nfeatures = 3000)
  scRNA.anchors <- FindIntegrationAnchors(object.list = scRNAlist, anchor.features = scRNA.features)
}
scRNA <- IntegrateData(anchorset = scRNA.anchors, dims = 1:30)
scRNA
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, features = rownames(scRNA))
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA, ndims = 50)
scRNA <- FindNeighbors(scRNA, dims = 1:50)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = 1:50)
scRNA <- RunTSNE(scRNA, dims = 1:50,check_duplicates = FALSE)
head(Idents(scRNA), 5)
save(scRNA,file="PBMC_scRNA.RData")

