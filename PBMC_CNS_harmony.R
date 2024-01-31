################################################################################
###-----harmony_integration----
################################################################################
library(Seurat);library(harmony);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

scRNA_PBMC = PBMC_scRNA
snRNA_EC = EC_snRNA
snRNA_HIP = HIP_snRNA
snRNA_PFC = PFC_snRNA
snRNA_MTG = MTG_snRNA
Idents(scRNA)="SampleID"
table(Idents(scRNA))
snRNA_PFC = subset(scRNA,downsample=10000)
rm(scRNA)

scRNA_PBMC[["Region"]] <- "PBMC"
scRNA_PBMC$celltype_Region_Group_Raw <- paste(scRNA_PBMC$celltype, scRNA_PBMC$Region, scRNA_PBMC$Diagnosis, sep = "_")
snRNA_EC$celltype_Region_Group_Raw <- paste(snRNA_EC$celltype, snRNA_EC$Region, snRNA_EC$Diagnosis, sep = "_")
snRNA_HIP$celltype_Region_Group_Raw <- paste(snRNA_HIP$celltype, snRNA_HIP$Region, snRNA_HIP$Diagnosis, sep = "_")
snRNA_PFC$celltype_Region_Group_Raw <- paste(snRNA_PFC$celltype, snRNA_PFC$Region, snRNA_PFC$Diagnosis, sep = "_")
snRNA_MTG$celltype_Region_Group_Raw <- paste(snRNA_MTG$celltype, snRNA_MTG$Region, snRNA_MTG$Diagnosis, sep = "_")

scRNA_snRNA <- merge(snRNA_EC, y=scRNA_PBMC) ##merge
scRNA_snRNA <- merge(snRNA_HIP, y=scRNA_PBMC) ##merge
scRNA_snRNA <- merge(snRNA_PFC, y=scRNA_PBMC) ##merge
scRNA_snRNA <- merge(snRNA_MTG, y=scRNA_PBMC) ##merge
rm(scRNA_PBMC,snRNA_PFC,snRNA_MTG)

table(scRNA_snRNA$orig.ident)
table(scRNA_snRNA$SampleID)
colnames(scRNA_snRNA@meta.data) #
cellinfo <- subset(scRNA_snRNA@meta.data, select= c("SampleID","Age","Sex","Region","Diagnosis","celltype","percent.mt","celltype_Region_Group_Raw")) #提取感兴趣meta信息
scRNA_snRNA <- CreateSeuratObject(scRNA_snRNA@assays$RNA@counts, meta.data = cellinfo)
###
scRNA_snRNA <- SCTransform(scRNA_snRNA)
###-
scRNA_snRNA <- RunPCA(scRNA_snRNA, npcs=50, verbose=FALSE)
###
scRNA_snRNA <- RunHarmony(scRNA_snRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
ElbowPlot(scRNA_snRNA, ndims = 20)
scRNA_snRNA <- RunTSNE(scRNA_snRNA, reduction="harmony",dims=1:30) %>% RunUMAP(reduction="harmony", dims = 1:30)
DimPlot(scRNA_snRNA,group.by = "SampleID",reduction = "umap")
DimPlot(scRNA_snRNA,group.by = "celltype",reduction = "umap",label = TRUE)
DimPlot(scRNA_snRNA,group.by = "SampleID", reduction = "umap", split.by = "Region")

DefaultAssay(scRNA_snRNA) <- "SCT"
scRNA_snRNA <- FindNeighbors(scRNA_snRNA, dims = 1:30) %>% FindClusters(dims = 1:30, resolution = 1.2)
#load("ref_Hematopoietic.RData")
#scRNA_snRNA <- cell_identify(scRNA_snRNA, ref_Hematopoietic)
head(Idents(scRNA_snRNA), 5)
save(scRNA_snRNA,file="PBMC_PFC4.RData")
save(scRNA_snRNA,file="PBMC_PFC2.RData")
save(scRNA_snRNA,file="PBMC_PFC1.RData")
save(scRNA_snRNA,file="PBMC_EC.RData")
save(scRNA_snRNA,file="PBMC_HIP.RData")

table(scRNA_snRNA@meta.data[["seurat_clusters"]], scRNA_snRNA@meta.data[["Sex"]])
table(scRNA_snRNA@meta.data[["seurat_clusters"]], scRNA_snRNA@meta.data[["celltype"]])
table(scRNA_snRNA@meta.data[["Age"]],scRNA_snRNA@meta.data[["seurat_clusters"]])
table(scRNA_snRNA@meta.data[["seurat_clusters"]],scRNA_snRNA@meta.data[["Age"]])

{#tsne
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "SampleID", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "SampleID", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "SampleID", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_SampleID_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Region", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Region_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Region", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Region_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "Region", groupFacet = 'Region', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Region_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Diagnosis", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Diagnosis_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Diagnosis", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Diagnosis_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "Diagnosis", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Diagnosis_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Age", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Age_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Age", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Age_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "Age", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Age_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_cluster_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  #umap
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "SampleID", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "SampleID", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "SampleID", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_SampleID_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Diagnosis", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Diagnosis_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Diagnosis", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Diagnosis_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "Diagnosis", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Diagnosis_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Age", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Age_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Age", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Age_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "Age", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Age_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_cluster_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  rm('p1')
}
#
{#tsne
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "SampleID", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "SampleID", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "SampleID", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_SampleID_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Diagnosis", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Diagnosis_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Diagnosis", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Diagnosis_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "Diagnosis", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Diagnosis_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Age", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Age_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "Age", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Age_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "Age", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Age_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = T,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_1.pdf", plot = p1, device = 'pdf', width = 14.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_2.pdf", plot = p1, device = 'pdf', width = 15, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'tsne',clusterCol = "celltype", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  #umap
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "SampleID", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "SampleID", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "SampleID", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_SampleID_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Diagnosis", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Diagnosis_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Diagnosis", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Diagnosis_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "Diagnosis", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Diagnosis_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Age", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Age_1.pdf", plot = p1, device = 'pdf', width = 11.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "Age", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Age_2.pdf", plot = p1, device = 'pdf', width = 12, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "Age", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Age_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 14.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA_snRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_2.pdf", plot = p1, device = 'pdf', width = 15, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA_snRNA,reduction = 'umap',clusterCol = "celltype", groupFacet = 'Diagnosis', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  rm('p1')
}

###Dotplot
DefaultAssay(scRNA_snRNA) <- "RNA"
scRNA_snRNA@meta.data[["celltype"]]<-factor(scRNA_snRNA@meta.data[["celltype"]], levels=c("RGL","IPC","Microglia","OPC","NFOL","MOL",
                                                                              "ImN","DG_ExN","nonDG_ExN","InN","CR",
                                                                              "End","Unk"))
scRNA@meta.data[["celltype"]]<-factor(scRNA@meta.data[["celltype"]], levels=c("Astro","B_Memory","B_Naive","B_Plasma",
                                                                              "CCR6_Th","CCR7_Th","CMP","DC","EPC","ExN",
                                                                              "FOXP3negCTLA4pos_Treg","FOXP3posCTLA4pos_Treg",
                                                                              "GZMB_KLRB1_Tc","GZMB_Tc","GZMK_Tc","KLRB1_Tc",
                                                                              "Micro","Monocyte","NK","Oligo","OPC","Peri","Platelets","T_cell"))

table(scRNA_snRNA@meta.data$celltype)
markers.to.plot <- c("CD19","MS4A1","AIM2","IGHD","JCHAIN","IGKC","CD4","IL7R",	"CCR6",	"CCR7",
                     "CD8A","CD3D","GZMB","GZMK","KLRB1","FOXP3",	"CTLA4","TIGIT","GNLY","KLRF1",
                     "CD3E","CD3G","PPBP","PF4","LYZ","S100A9","SOX4","ZNF385D","HBA1","HBA2","CLNK","WDFY4",
                     "GFAP","AQP4","SOX2","SLC1A3","CNTNAP2","SYT1","RBFOX1","SNAP25","SLC17A7","MAP1B","SATB2","CAMK2A","STMN2","GAD1","GAD2","DCX",
                     "CSF1R","P2RY12","CD74","CX3CR1","MOBP","PLP1","MBP","BCAS1","OPALIN","VCAN","SOX6","PDGFRA",
                     "FLT1","CLDN5","DCN","COL1A2","PTPRC","SKAP1","CD247")
pdf("cluster_marker_1.pdf", width = 20,height = 10)
DotPlot(scRNA_snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("cluster_marker_2.pdf", width = 20,height = 10)
DotPlot(scRNA_snRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
pdf("celltype_marker_1.pdf", width = 20,height = 7)
DotPlot(scRNA_snRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("celltype_marker_2.pdf", width = 20,height = 7)
DotPlot(scRNA_snRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

##
table(scRNA_snRNA$stage)  
av<-AverageExpression(scRNA_snRNA,group.by = "stage", assays = "SCT")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_stage_1.pdf", width = 5.5,height = 5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_stage_1.csv")
#
table(scRNA_snRNA$seurat_clusters)  
av<-AverageExpression(scRNA_snRNA,group.by = "seurat_clusters", assays = "SCT")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_seurat_clusters_1.pdf", width = 5.5,height = 5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_seurat_clusters_1.csv")
#
table(scRNA_snRNA$celltype)  
av<-AverageExpression(scRNA_snRNA,group.by = "celltype", assays = "SCT")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype.pdf", width = 10,height = 9)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv") 
##
save(scRNA_snRNA,file="scRNA_snRNA.RData")






