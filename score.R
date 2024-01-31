################################################################################
###---score
################################################################################
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr)
DefaultAssay(scRNA)="RNA"
blue_gene<-read_xlsx("gene.xlsx")
gene<-as.list(blue_gene)
AB=scRNA
AB<-AddModuleScore(AB, features = gene, ctrl = 100, name = "gene")
colnames(AB@meta.data)[18]<-"blue_Score"
##
P1<-VlnPlot(AB, features = 'blue_Score',pt.size = 0.01)
ggsave(filename = "blue_1.pdf", plot = P1, device = 'pdf', width = 20, height = 15, units = 'cm')
rm('P1')  
###
P1<-RidgePlot(AB, features = 'blue_Score', ncol = 1) 
ggsave(filename = "blue_2.pdf", plot = P1, device = 'pdf', width = 18, height = 15, units = 'cm')
rm('P1')  
##
data<-FetchData(AB, vars = c("celltype","blue_Score"))
P1<-ggplot(data,aes(celltype, blue_Score))+
  geom_boxplot()+theme_bw()+RotatedAxis()
ggsave(filename = "blue_3.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 

##Distribution
b<-FetchData(AB, vars = c("Diagnosis","blue_Score"))
b[["Diagnosis"]]<-factor(b[["Diagnosis"]], levels=c("Control","Early_AD","Late_AD"))
P1<-ggdensity(b, x = "blue_Score",
          add = "mean", rug = TRUE,
          color = "Diagnosis", fill = "Diagnosis",
          palette = c("#00AFBB", "#E7B800","#8B658B"))
ggsave(filename = "blue_4.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 
P1<-gghistogram(b, x = "blue_Score",bins = 30, add = "mean", rug = TRUE,
            color = "Diagnosis", fill = "Diagnosis",
            palette = c("#00AFBB", "#E7B800","#8B658B"))
ggsave(filename = "blue_5.pdf", plot = P1, device = 'pdf', width = 15, height = 15, units = 'cm')
rm('P1') 

##
b<-FetchData(AB, vars = c("Diagnosis","blue_Score"))
my_comparisons <- list( c("Control", "Early_AD"), c("Control", "Late_AD"), c("Early_AD", "Late_AD"))
ggboxplot(b, x = "Diagnosis", y = "blue_Score",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "Diagnosis", palette = "npg")+stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=45,size=6))
ggsave("blue_6.pdf",width = 15,height = 15,units = "cm")

ggviolin(b, x = "Diagnosis", y = "blue_Score", fill = "Diagnosis",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)
ggsave("blue_7.pdf",width = 15,height = 15,units = "cm")

##palette = c("#00AFBB", "#E7B800", "#FC4E07")

##
b<-FetchData(AB, vars = c("Diagnosis","blue_Score","celltype"))
b[["Diagnosis"]]<-factor(b[["Diagnosis"]], levels=c("Control","Early_AD","Late_AD"))
ggdensity(b, x="blue_Score", facet.by = "celltype",y="..density..", combine = TRUE,
          xlab = "blue_Score", add = "median", rug = TRUE, color = "Diagnosis", 
          fill = "Diagnosis", palette = "jco",ncol=4)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=45,size=6))
ggsave("blue_8.pdf",width = 15,height = 15,units = "cm")
##
ggboxplot(b, x = "Diagnosis", y = "blue_Score",facet.by = "celltype",combine = TRUE,add = "jitter", 
          add.params = list(size=0.1, jitter=0.2),label.select = list(top.up=2, top.down=2), 
          font.label = list(size=9, face="italic"), repel = TRUE,
          color = "Diagnosis", palette = "npg",ncol=4)+stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1, angle=45,size=8))
ggsave("blue_9.pdf",width = 20,height = 25,units = "cm")
##
ggviolin(b, x = "Diagnosis", y = "blue_Score",facet.by = "celltype", fill = "Diagnosis",
         palette = "jitter",
         add = "boxplot", add.params = list(fill = "white"),size = 0)+
  stat_compare_means(comparisons = my_comparisons)+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1, angle=45,size=8))
ggsave("blue_10.pdf",width = 20,height = 20,units = "cm")


