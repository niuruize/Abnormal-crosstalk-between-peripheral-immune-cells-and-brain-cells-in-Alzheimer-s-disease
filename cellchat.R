#loading R packages
library(Seurat);library(tidyverse);library(cowplot);library(patchwork)
library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db)
library(MySeuratWrappers);library(scRNAtoolVis);library(readxl)

################################################################################
###---cellchat
################################################################################
library(CellChat); library(ggplot2); library(ggalluvial); library(Seurat)
##
Control <- scRNA[,scRNA@meta.data$Diagnosis %in% c("Control")]
Early_AD <- scRNA[,scRNA@meta.data$Diagnosis %in% c("Early_AD")]
Late_AD <- scRNA[,scRNA@meta.data$Diagnosis %in% c("Late_AD")]

{data.input  <- Control@assays[["RNA"]]@data
  celltype  <- Control@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Control$celltype, row.names = names(Control$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Control$celltype, row.names = names(identify))}
##
{data.input  <- Early_AD@assays[["RNA"]]@data
  celltype  <- Early_AD@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Early_AD$celltype, row.names = names(Early_AD$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Early_AD$celltype, row.names = names(identify))
}
##
{data.input  <- Late_AD@assays[["RNA"]]@data
  celltype  <- Late_AD@meta.data[["celltype"]]
  data.input[1:4,1:4]
  identity = data.frame(group = Late_AD$celltype, row.names = names(Late_AD$celltype)) # create a dataframe consisting of the cell labels
  head(identity)
  unique(identity$group) # check the cell labels
  table(identity$group)
  meta <- data.frame(labels = Late_AD$celltype, row.names = names(identify))}
##
cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels")
cellchat
summary(cellchat)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
head(cellchat@meta)
## set "labels" as default cell identity
cellchat <- setIdent(cellchat, ident.use = "labels") 
# show factor levels of the cell labels
levels(cellchat@idents) 
groupSize <- as.numeric(table(cellchat@idents))
table(cellchat@idents)
##
CellChatDB <- CellChatDB.human 
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
# set the used database in the object
#unique(CellChatDB$interaction$annotation)
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
#CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") 
#CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") 
cellchat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) 
######plan("multiprocess", workers = 4) # do parallel  
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat <- computeCommunProb(cellchat,raw.use = T, population.size = T,type = "truncatedMean",trim = 0.1)  
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat,slot.name = "netP")
#cellchat <- computeNetSimilarity(cellchat, type = "functional")
#cellchat <- netEmbedding(cellchat, type = "functional")
#cellchat <- netClustering(cellchat, type = "functional")
#cellchat <- computeNetSimilarity(cellchat, type = "structural")
#cellchat <- netEmbedding(cellchat, type = "structural")
#cellchat <- netClustering(cellchat, type = "structural")

###------------------
Control <- cellchat
saveRDS(Control, "Control.rds")
Early_AD <- cellchat
saveRDS(Early_AD, "Early_AD.rds")
Late_AD <- cellchat
saveRDS(Late_AD, "Late_AD.rds")

cellchat_all<-data.frame(Control=sort(Control@netP$pathways),Early_AD=sort(Early_AD@netP$pathways),Late_AD=sort(Late_AD@netP$pathways))
cellchat_Control=data_frame(Control=sort(Control@netP$pathways))
cellchat_Early_AD=data_frame(Early_AD=sort(Early_AD@netP$pathways))
cellchat_Late_AD=data_frame(Late_AD=sort(Late_AD@netP$pathways))

cellchat_all <-dplyr::bind_rows(cellchat_Control,cellchat_Early_AD,cellchat_Late_AD)
write.csv(cellchat_all,"cellchat_all.csv") 
cellchat_all1<-data.frame(cellchat_all=intersect(x=Control@netP$pathways, y = c(Early_AD@netP$pathways, Late_AD@netP$pathways)))
write.csv(cellchat_all1,"cellchat_all1.csv") 

###------------------
## merge
con.list <- list(Control=Control,Early_AD=Early_AD,Late_AD=Late_AD)
cellchat <- mergeCellChat(con.list,add.names = names(con.list),cell.prefix = T)
saveRDS(cellchat, "cellchat.rds")

# count-weight
gg1 <- compareInteractions(cellchat, show.legend = F, color.use = c("#00AFBB", "#E7B800","#8B658B") ,group = c(1,2,3),measure = "count")+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1,angle=45,size=8))
gg2 <- compareInteractions(cellchat, show.legend = F, color.use = c("#00AFBB", "#E7B800","#8B658B") ,group = c(1,2,3),measure = "weight")+
  theme(axis.text.x=element_text(vjust = 1, hjust = 1,angle=45,size=8))
ggsave("overview_number_strength.pdf", gg1+gg2, width = 3,height = 2)











