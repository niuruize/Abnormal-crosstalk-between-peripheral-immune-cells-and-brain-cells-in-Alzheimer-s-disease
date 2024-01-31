################################################################################
###---HIP
################################################################################
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

table(scRNA$orig.ident)
table(scRNA$SampleID)
table(scRNA$group)
colnames(scRNA@meta.data)
scRNAlist <- SplitObject(scRNA, split.by = "SampleID")
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst",nfeatures = 3000,verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = scRNAlist,nfeatures = 3000)
scRNAlist <- lapply(X = scRNAlist, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
##
anchors <- FindIntegrationAnchors(object.list = scRNAlist, reference = c(5,10,15,8), reduction = "rpca", dims = 1:50)
scRNA.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
scRNA=scRNA.integrated
save(scRNA,file="PFC1_snRNA.RData")
##first
DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30, check_duplicates = FALSE) 
head(Idents(scRNA), 5)
##second
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5","6","7","8","9","10",
                                            "11","12","14","16","17","18","19","20",
                                            "21","22","23","24","25","26","32","33","34")]

DefaultAssay(scRNA) <- "integrated"
# Run the standard workflow for visualization and clustering
scRNA <- ScaleData(scRNA, verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)
ElbowPlot(scRNA)
scRNA <- FindNeighbors(scRNA, dims = 1:30)
scRNA <- FindClusters(scRNA, resolution = 0.4)
scRNA <- RunUMAP(scRNA, dims = 1:30)
scRNA <- RunTSNE(scRNA, dims = 1:30, check_duplicates = FALSE) 
head(Idents(scRNA), 5)

{#tsne
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "datasets", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_datasets_1.pdf", plot = p1, device = 'pdf', width = 12.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "datasets", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_datasets_2.pdf", plot = p1, device = 'pdf', width = 13, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'datasets', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_datasets_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "stage", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_stage_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "stage", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_stage_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'stage', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_stage_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sex", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Sex_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sex", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Sex_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'Sex', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Sex_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_1.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_2.pdf", plot = p1, device = 'pdf', width = 12.5, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_cluster_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_1.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_2.pdf", plot = p1, device = 'pdf', width = 12.5, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype_3.pdf", plot = p1, device = 'pdf', width = 20, height = 10, units = 'cm')
  
  #umap
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "datasets", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_datasets_1.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "datasets", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_datasets_2.pdf", plot = p1, device = 'pdf', width = 12.5, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'datasets', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_datasets_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "stage", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_stage_1.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "stage", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_stage_2.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'stage', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_stage_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sex", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Sex_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sex", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Sex_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'Sex', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Sex_3.pdf", plot = p1, device = 'pdf', width = 40, height = 10, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_cluster_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  rm('p1')
}
#
Idents(scRNA)="seurat_clusters"
{#tsne
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sex", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sex", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "Sex", groupFacet = 'Sex', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_SampleID_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "group", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Group_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "group", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Group_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "group", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Group_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_1.pdf", plot = p1, device = 'pdf', width = 30, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_2.pdf", plot = p1, device = 'pdf', width = 30, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_cluster_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  #umap
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sex", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sex", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "Sex", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_SampleID_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "group", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Group_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "group", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Group_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "group", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Group_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_1.pdf", plot = p1, device = 'pdf', width = 28, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_2.pdf", plot = p1, device = 'pdf', width = 28, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_cluster_3.pdf", plot = p1, device = 'pdf', width = 20, height = 9, units = 'cm')
  
  rm('p1')
}
#
Idents(scRNA)="celltype"
{#tsne
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = T,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_1.pdf", plot = p1, device = 'pdf', width = 28, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_2.pdf", plot = p1, device = 'pdf', width = 28, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype_3.pdf", plot = p1, device = 'pdf', width = 30, height = 18, units = 'cm')
  
  #umap
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 30, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_2.pdf", plot = p1, device = 'pdf', width = 30, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype", groupFacet = 'group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype_3.pdf", plot = p1, device = 'pdf', width = 30, height = 18, units = 'cm')
  
  rm('p1')
}

###Dotplot
DefaultAssay(scRNA) <- "RNA"
markers.to.plot <- c("MOBP","PLP1","MBP","OPALIN","CNTNAP2","SYT1",
                     "RBFOX1","SNAP25","SLC17A7","MAP1B","SATB2","CAMK2A","STMN2","PROX1","GAD1","GAD2","DCX",
                     "FLT1","CLDN5","DCN","COL1A2","VCAN","SOX6","PDGFRA","GFAP","AQP4","SOX2","SLC1A3","GPC5","SLC1A2","ATP1B2",
                     "CSF1R","C3","P2RY12","CD74","CX3CR1","PTPRC","LPAR6","DOCK8","CCK","LAMP5","SV2C","CNR1","SST","VIP","NOS1","PVALB","CLIC6")
markers.to.plot <- c("GFAP","AQP4","SOX2","SLC1A3","NETO1","PROX1","BCL11B","DCX","SEMA3C","MPPED1","CPNE4","SATB2","PCP4","RGS4",
                     "NNAT","STMN2","CAMK2A","GAD1","GAD2","MOBP","MOG","PDGFRA","OLIG2","CX3CR1","PTPRC","CCK","LAMP5","SV2C","CNR1",
                     "SST","CALB2","PVALB","RELN","DNAH9","FOXJ1","CFAP54","FLI1","FN1","DCN","COL1A2")
markers.to.plot <- c("RBFOX1","SNAP25","CAMK2A","SATB2","GAD1","GAD2","GFAP","AQP4","VCAN","SOX6","PLP1","MBP",
                     "CSF1R","P2RY12","CD74","CX3CR1","FLT1","CLDN5","DCN","COL1A2","PTPRC","CD247")
pdf("cluster_marker_1.pdf", width = 14,height = 12)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("cluster_marker_2.pdf", width = 14,height = 12)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
pdf("cluster_marker_3.pdf", width = 12,height = 12)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("cluster_marker_4.pdf", width = 12,height = 12)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
pdf("celltype_marker_1.pdf", width = 7.5,height = 3.5)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("celltype_marker_2.pdf", width = 7.5,height = 3.5)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

##
cluster_celltype <-  data.frame(readxl::read_xlsx("cluster_celltype.xlsx"))
#scRNA = scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")]
current.cluster.ids <- cluster_celltype$cluster
new.cluster.ids <- cluster_celltype$celltype1
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$seurat_clusters,scRNA$celltype)
Idents(scRNA) <- factor(Idents(scRNA), levels=c("Astrocyte","ExN","InN","OPC","Oligodendrocyte","Microglia","Endotheliocyte"))
scRNA@meta.data[["celltype"]]<-factor(scRNA@meta.data[["celltype"]], levels=c("ExN","InN","Astro","OPC","Oligo","Micro","Endo","Peri"))
scRNA$Diagnosis<-factor(scRNA$Diagnosis, levels=c("Control","AD"))

##
table(scRNA$seurat_clusters)  
av<-AverageExpression(scRNA,group.by = "seurat_clusters", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_seurat_clusters_1.pdf", width = 7.5,height = 7)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_seurat_clusters_1.csv")
#
table(scRNA$celltype)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype.pdf", width = 3,height = 2.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv")
#
table(scRNA$celltype2)  
av<-AverageExpression(scRNA,group.by = "celltype2", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype2.pdf", width = 7,height = 6.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype2.csv")
#
table(scRNA$celltype)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype.pdf", width = 8,height = 7.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv") 

##
Idents(scRNA)="seurat_clusters"
Idents(scRNA)="celltype"
prop.table(table(Idents(scRNA)))
table(Idents(scRNA), scRNA$Diagnosis)
Cellratio <- prop.table(table(Idents(scRNA), scRNA$Diagnosis), margin = 2)
Cellratio <- as.data.frame(Cellratio)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = my36colors)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(filename = "Cellratio_all.pdf", device = 'pdf', width = 10, height = 12, units = 'cm')

##
save(scRNA,file="PFC1_snRNA.RData")

