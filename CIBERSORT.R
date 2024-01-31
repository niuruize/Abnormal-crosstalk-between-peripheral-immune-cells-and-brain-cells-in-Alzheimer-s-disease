#loading R packages
library(CIBERSORT) #devtools::install_github("Moonerss/CIBERSORT")
#library(Mfuzz) #BiocManager::install('Mfuzz')
library(Seurat);library(progeny);library(tidyr);library(tibble);library(dplyr);library(viridisLite);library(ggplot2);library(ggpubr);library(ggsci);library(readxl)
################################################################################
###---CIBERSORT
################################################################################
# The input data were obtained from cross "PBMC_integration.R"
load("scRNA.RData")
age.averages <- AverageExpression(scRNA, group.by = "celltype")
head(age.averages[["integrated"]])
sig_matrix <- as.matrix(age.averages[["RNA"]])

## bulk expression matrix
exp_data <- data.frame(read_xlsx("exp_data_1.xlsx"))
for (i in 1:ncol(exp_data)){exp_data[,i][is.na(exp_data[,i])] <- 0}
colnames(exp_data)
expr_read = data.matrix(exp_data[,2:88]) 
rownames(expr_read) = exp_data$Symbol
expr_read[1:4,1:4]
##
results <- cibersort(sig_matrix, expr_read, perm = 1000,QN = T)
head(results[,1:4],n=20)
write.csv(results, 'TCGA_CIBERSORT.csv')
## metadata
group1 <- data.frame(sample=TCGAquery_SampleTypes(colnames(expData), typesample = c("NT")),group="NT")
group2 <- data.frame(sample=TCGAquery_SampleTypes(colnames(expData), typesample = c("TP")),group="TP")
use.clin <- rbind(group1,group2)
results1 <- data.frame(results)
results2 <- mutate(results1,sample = rownames(results1))
results3 <- merge(use.clin,  results2, by = "sample", all=TRUE)
rownames(results3) <- results3$sample
results3 <- results3[,-1]
res <- data.frame(results3[,1:8])%>%
  rownames_to_column("sample")%>%
  pivot_longer(cols = colnames(.)[3:9],
               names_to = "cell.type",
               values_to = 'value')
head(res,n=6)

ggplot(res,aes(cell.type,value,fill = group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1, size = 10,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 10,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(aes(group = group),label = "p.format",size=3,method = "kruskal.test")
ggsave(filename = "TCGA_细胞丰度_1.pdf",  device = 'pdf', width = 15, height = 10, units = 'cm')

normalize <- function(x) {
  if((max(x) - min(x)) == 0){
    return(mean(x))
  }else{
    return((x - min(x)) / (max(x) - min(x)))
  }
}

heat_map_res <- apply(results[,1:7], 2, normalize)
heat_map_res <- t(as.data.frame(heat_map_res))
mycol<-colorRampPalette(c("navy", "white", "firebrick3"))(100)
library(pheatmap)
p1=pheatmap(heat_map_res,color = mycol,cluster_rows = T,cluster_cols = T,cellwidth = 0.5, cellheight = 10,treeheight_row =10, treeheight_col=10,show_colnames = F,scale = 'none')
ggsave(filename = "TCGA_细胞丰度_2.pdf", plot = p1, device = 'pdf', width = 12, height = 5, units = 'cm')












