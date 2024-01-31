#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
###---DEGs
################################################################################
#---------------------------DEGs---------------------------###
library(ggplot2); library(Seurat); library(patchwork); library(cowplot)
colnames(scRNA@meta.data)
Idents(scRNA) <- "age"
scRNA$celltype_Diagnosis <- paste(scRNA$celltype, scRNA$Diagnosis, sep = "_")
table(scRNA$celltype_Diagnosis)
Idents(scRNA) <- "celltype_Diagnosis"
DefaultAssay(scRNA)="RNA"
###1
{
age.response <- FindMarkers(scRNA, ident.1 = c("B_Memory_Early_AD"), ident.2 = ("B_Memory_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Memory_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("B_Memory_Late_AD"), ident.2 = ("B_Memory_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Memory_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("B_Memory_Late_AD"), ident.2 = ("B_Memory_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Memory_Late_Early.csv")
}
###2
{
age.response <- FindMarkers(scRNA, ident.1 = c("B_Naive_Early_AD"), ident.2 = ("B_Naive_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Naive_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("B_Naive_Late_AD"), ident.2 = ("B_Naive_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Naive_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("B_Naive_Late_AD"), ident.2 = ("B_Naive_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Naive_Late_Early.csv")
}
###3
{
age.response <- FindMarkers(scRNA, ident.1 = c("B_Plasma_Early_AD"), ident.2 = ("B_Plasma_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Plasma_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("B_Plasma_Late_AD"), ident.2 = ("B_Plasma_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Plasma_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("B_Plasma_Late_AD"), ident.2 = ("B_Plasma_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="B_Plasma_Late_Early.csv")
}
###4
{
age.response <- FindMarkers(scRNA, ident.1 = c("CCR6_Th_Early_AD"), ident.2 = ("CCR6_Th_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CCR6_Th_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("CCR6_Th_Late_AD"), ident.2 = ("CCR6_Th_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CCR6_Th_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("CCR6_Th_Late_AD"), ident.2 = ("CCR6_Th_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CCR6_Th_Late_Early.csv")
}
###5
{
age.response <- FindMarkers(scRNA, ident.1 = c("CCR7_Th_Early_AD"), ident.2 = ("CCR7_Th_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CCR7_Th_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("CCR7_Th_Late_AD"), ident.2 = ("CCR7_Th_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CCR7_Th_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("CCR7_Th_Late_AD"), ident.2 = ("CCR7_Th_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CCR7_Th_Late_Early.csv")
}
###6
{
age.response <- FindMarkers(scRNA, ident.1 = c("FOXP3-CTLA4+_Treg_Early_AD"), ident.2 = ("FOXP3-CTLA4+_Treg_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="FOXP3-CTLA4+_Treg_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("FOXP3-CTLA4+_Treg_Late_AD"), ident.2 = ("FOXP3-CTLA4+_Treg_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="FOXP3-CTLA4+_Treg_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("FOXP3-CTLA4+_Treg_Late_AD"), ident.2 = ("FOXP3-CTLA4+_Treg_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="FOXP3-CTLA4+_Treg_Late_Early.csv")
}
###7
{
age.response <- FindMarkers(scRNA, ident.1 = c("FOXP3+CTLA4+_Treg_Early_AD"), ident.2 = ("FOXP3+CTLA4+_Treg_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="FOXP3+CTLA4+_Treg_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("FOXP3+CTLA4+_Treg_Late_AD"), ident.2 = ("FOXP3+CTLA4+_Treg_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="FOXP3+CTLA4+_Treg_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("FOXP3+CTLA4+_Treg_Late_AD"), ident.2 = ("FOXP3+CTLA4+_Treg_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="FOXP3+CTLA4+_Treg_Late_Early.csv")
}
###8
{
age.response <- FindMarkers(scRNA, ident.1 = c("GZMB_KLRB1_Tc_Early_AD"), ident.2 = ("GZMB_KLRB1_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMB_KLRB1_Tc_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("GZMB_KLRB1_Tc_Late_AD"), ident.2 = ("GZMB_KLRB1_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMB_KLRB1_Tc_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("GZMB_KLRB1_Tc_Late_AD"), ident.2 = ("GZMB_KLRB1_Tc_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMB_KLRB1_Tc_Late_Early.csv")
}
###9
{
age.response <- FindMarkers(scRNA, ident.1 = c("GZMB_Tc_Early_AD"), ident.2 = ("GZMB_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMB_Tc_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("GZMB_Tc_Late_AD"), ident.2 = ("GZMB_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMB_Tc_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("GZMB_Tc_Late_AD"), ident.2 = ("GZMB_Tc_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMB_Tc_Late_Early.csv")
}
###10
{
age.response <- FindMarkers(scRNA, ident.1 = c("GZMK_Tc_Early_AD"), ident.2 = ("GZMK_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMK_Tc_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("GZMK_Tc_Late_AD"), ident.2 = ("GZMK_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMK_Tc_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("GZMK_Tc_Late_AD"), ident.2 = ("GZMK_Tc_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="GZMK_Tc_Late_Early.csv")
}
###11
{
age.response <- FindMarkers(scRNA, ident.1 = c("KLRB1_Tc_Early_AD"), ident.2 = ("KLRB1_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="KLRB1_Tc_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("KLRB1_Tc_Late_AD"), ident.2 = ("KLRB1_Tc_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="KLRB1_Tc_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("KLRB1_Tc_Late_AD"), ident.2 = ("KLRB1_Tc_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="KLRB1_Tc_Late_Early.csv")
}
###12
{
age.response <- FindMarkers(scRNA, ident.1 = c("Platelets_Early_AD"), ident.2 = ("Platelets_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="Platelets_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("Platelets_Late_AD"), ident.2 = ("Platelets_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="Platelets_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("Platelets_Late_AD"), ident.2 = ("Platelets_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="Platelets_Late_Early.csv")
}
###13
{
age.response <- FindMarkers(scRNA, ident.1 = c("NK_Early_AD"), ident.2 = ("NK_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="NK_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("NK_Late_AD"), ident.2 = ("NK_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="NK_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("NK_Late_AD"), ident.2 = ("NK_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="NK_Late_Early.csv")
}
###14
{
age.response <- FindMarkers(scRNA, ident.1 = c("Monocyte_Early_AD"), ident.2 = ("Monocyte_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="Monocyte_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("Monocyte_Late_AD"), ident.2 = ("Monocyte_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="Monocyte_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("Monocyte_Late_AD"), ident.2 = ("Monocyte_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="Monocyte_Late_Early.csv")
}
###15
{
age.response <- FindMarkers(scRNA, ident.1 = c("CMP_Early_AD"), ident.2 = ("CMP_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CMP_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("CMP_Late_AD"), ident.2 = ("CMP_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CMP_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("CMP_Late_AD"), ident.2 = ("CMP_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="CMP_Late_Early.csv")
}
###16
{
age.response <- FindMarkers(scRNA, ident.1 = c("EPC_Early_AD"), ident.2 = ("EPC_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="EPC_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("EPC_Late_AD"), ident.2 = ("EPC_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="EPC_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("EPC_Late_AD"), ident.2 = ("EPC_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="EPC_Late_Early.csv")
}
###17
{
age.response <- FindMarkers(scRNA, ident.1 = c("DC_Early_AD"), ident.2 = ("DC_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="DC_Early_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("DC_Late_AD"), ident.2 = ("DC_Control"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="DC_Late_Control.csv")
age.response <- FindMarkers(scRNA, ident.1 = c("DC_Late_AD"), ident.2 = ("DC_Early_AD"),verbose = FALSE)
head(age.response, n = 10)
write.csv(age.response,file="DC_Late_Early.csv")
}











