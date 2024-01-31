################################################################################
###---NicheNet
################################################################################
#devtools::install_github("saeyslab/nichenetr")
library(nichenetr);library(RColorBrewer);library(tidyverse);library(Seurat);library(circlize)
##
ligand_target_matrix = readRDS("/Users/niuruize/Downloads/scRNA/NicheNet/ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5]
lr_network = readRDS("/Users/niuruize/Downloads/scRNA/NicheNet/lr_network.rds")
head(lr_network)
weighted_networks = readRDS("/Users/niuruize/Downloads/scRNA/NicheNet/weighted_networks.rds")
head(weighted_networks$lr_sig)
head(weighted_networks$gr)
# Seurat
Idents(scRNA)="celltype"
DefaultAssay(scRNA)="RNA"
seuratObj = scRNA
seuratObj@meta.data %>% head()
seuratObj@meta.data$celltype %>% table()
DimPlot(seuratObj, reduction = "tsne")
seuratObj@meta.data$Diagnosis %>% table()
DimPlot(seuratObj, reduction = "tsne", group.by = "Diagnosis")
###----NicheNet----####
colnames(seuratObj@meta.data)
table(seuratObj$Diagnosis)
table(seuratObj$celltype)
#"B_Memory","B_Naive","B_Plasma","CCR6_Th","CCR7_Th","GZMB_Tc","GZMK_Tc","KLRB1_Tc","GZMB_KLRB1_Tc","FOXP3posCTLA4pos_Treg","FOXP3negCTLA4pos_Treg",
#"NK","Platelets","Monocyte","CMP","EPC","DC"
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = c("GZMK_Tc"),
  condition_colname = "Diagnosis", condition_oi = c("Early_AD","Late_AD"), condition_reference = "Control", 
  sender = c("B_Memory","B_Naive","B_Plasma",
             "CCR6_Th","CCR7_Th","GZMB_Tc","GZMB_KLRB1_Tc","FOXP3posCTLA4pos_Treg","FOXP3negCTLA4pos_Treg",
             "NK","Platelets","Monocyte","CMP","EPC","DC","KLRB1_Tc"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "human", assay_oi = "RNA")
##
ligand_activities = nichenet_output$ligand_activities
write.csv(ligand_activities,"ligand_activities.csv") 
##
nichenet_output$top_ligands
nichenet_output$ligand_expression_dotplot
p = DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
ggsave("top20_ligands_1.pdf", p, width = 9, height = 4.5)
##
p = DotPlot(seuratObj, features = c("MIF", "CALM1", "ICAM2","HLA-A","PTPRC","GSTP1","ITGB1","CD99","YARS","PECAM1","F11R","TGFB1","CD226") %>% rev(), cols = "RdYlBu") + RotatedAxis()
ggsave("WGCNA_blue_ligands_1.pdf", p, width = 8, height = 4.5)
##
p = DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), split.by = "Diagnosis", cols = "RdYlBu") + RotatedAxis()
ggsave("top20_ligands_compare_1.pdf", p, width = 10, height = 12)
p = DotPlot(seuratObj, features = c("MIF", "CALM1", "ICAM2","HLA-A","PTPRC","GSTP1","ITGB1","CD99","YARS","PECAM1","F11R","TGFB1","CD226") %>% rev(), split.by = "Diagnosis", cols = "RdYlBu") + RotatedAxis()
ggsave("WGCNA_blue_ligands_compare_1.pdf", p, width = 8, height = 12)
##
#VlnPlot(seuratObj, features = nichenet_output$top_ligands, split.by = "Diagnosis", pt.size = 0.01, combine = T)
p = VlnPlot(seuratObj, features = c("MIF"), split.by = "Diagnosis", pt.size = 0.01, combine = T)
##
p = nichenet_output$ligand_target_heatmap
ggsave("Heatmap_ligand_target_1.pdf", p, width = 1.8, height = 3.5)
##
p = nichenet_output$ligand_target_heatmap + scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + xlab("KLRB1_Tc_GZMK_Tc") + ylab("Prioritized immmune cell ligands")
##
nichenet_output$ligand_differential_expression_heatmap
##
nichenet_output$top_targets
##
p = DotPlot(seuratObj %>% subset(idents = c("GZMK_Tc")), features = nichenet_output$top_targets %>% rev(), cols = "RdYlBu",split.by = "Diagnosis") + RotatedAxis()
ggsave("Target_expression_dotplot_1.pdf", p, width = 6, height = 1.9)
##
p = VlnPlot(seuratObj %>% subset(idents=c("GZMK_Tc")), features = c("FOS","JUN","JUNB","HBB","IL7R","CCL4L2"), ncol = 6,split.by = "Diagnosis", pt.size = 0.01, combine = T)
ggsave("Target_expression_vlnplot_1.pdf", p, width = 15, height = 4)
p = VlnPlot(seuratObj %>% subset(idents=c("B_Memory","B_Naive","B_Plasma","CCR6_Th","CCR7_Th","GZMB_Tc","GZMK_Tc","KLRB1_Tc","GZMB_KLRB1_Tc","FOXP3posCTLA4pos_Treg","FOXP3negCTLA4pos_Treg",
"NK","Platelets","Monocyte","CMP","EPC","DC")), features = c("FOS","JUN","JUNB","HBB","IL7R","CCL4L2"), ncol = 3,split.by = "Diagnosis", pt.size = 0.01, combine = T)
ggsave("Target_expression_vlnplot_2.pdf", p, width = 20, height = 8)
##
p = nichenet_output$ligand_activity_target_heatmap
ggsave("Heatmap_ligand_activity_target_1.pdf", p, width = 20, height = 8)
##
ligand_target_df=nichenet_output$ligand_target_df # weight column = regulatory potential
write.csv(ligand_target_df,"ligand_target_df.csv") 
ligand_target_matrix=nichenet_output$ligand_target_matrix
write.csv(ligand_target_matrix,"lligand_target_matrix.csv")
##
p = nichenet_output$ligand_receptor_heatmap
ggsave("Heatmap_ligand_receptor_1.pdf", p, width = 5.9, height = 4.8)
##
ligand_receptor_matrix = nichenet_output$ligand_receptor_matrix
write.csv(ligand_receptor_matrix,"ligand_receptor_matrix.csv")
ligand_receptor_df = nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction
write.csv(ligand_receptor_df,"ligand_receptor_df.csv") 
##
nichenet_output$top_receptors
p = DotPlot(seuratObj %>% subset(idents = c("GZMK_Tc")), features = nichenet_output$top_receptors, cols = "RdYlBu",split.by = "Diagnosis") + RotatedAxis()
ggsave("Receprots_expression_dotplot_1.pdf", p, width = 12, height = 2)
##
p = VlnPlot(seuratObj %>% subset(idents=c("GZMK_Tc")), features = c("TNFRSF14","CD74","ITGA4","RPSA","TNFRSF1B","ITGB1","TGFBR2","TGFBR3"), ncol = 8,split.by = "Diagnosis", pt.size = 0.01, combine = T)
ggsave("Receptors_expression_vlnplot_1.pdf", p, width = 20, height = 4)
##
p = nichenet_output$ligand_receptor_heatmap_bonafide
ggsave("Heatmap_ligand_receptor_bonafide_1.pdf", p, width = 2.1, height = 2.8)
##
ligand_receptor_matrix_bonafide=nichenet_output$ligand_receptor_matrix_bonafide
write.csv(ligand_receptor_matrix_bonafide,"ligand_receptor_matrix_bonafide.csv")
ligand_receptor_df_bonafide=nichenet_output$ligand_receptor_df_bonafide
write.csv(ligand_receptor_df_bonafide,"ligand_receptor_df_bonafide.csv")

###---circos
#avg_expression_ligands = AverageExpression(seuratObj %>% subset(subset = Diagnosis == "Early_AD"),features = nichenet_output$top_ligands)
avg_expression_ligands = AverageExpression(seuratObj,features = nichenet_output$top_ligands)
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
  }) %>% t()
sender_ligand_assignment[1:4,1:4]
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x)>0})
names(sender_ligand_assignment)
(sender_ligand_assignment)
all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = nichenet_output$top_ligands %>% setdiff(unique_ligands)

Monocyte_specific_ligands = sender_ligand_assignment$Monocyte %>% names() %>% setdiff(general_ligands)
B_Memory_specific_ligands = sender_ligand_assignment$B_Memory %>% names() %>% setdiff(general_ligands)
B_Plasma_specific_ligands = sender_ligand_assignment$B_Plasma %>% names() %>% setdiff(general_ligands)
CCR6_Th_specific_ligands = sender_ligand_assignment$CCR6_Th %>% names() %>% setdiff(general_ligands)
CCR7_Th_specific_ligands = sender_ligand_assignment$CCR7_Th %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("Monocyte", times = Monocyte_specific_ligands %>% length()),
                  rep("B_Memory", times = B_Memory_specific_ligands %>% length()),
                  rep("B_Plasma", times = B_Plasma_specific_ligands %>% length()),
                  rep("CCR6_Th",  times = CCR6_Th_specific_ligands %>% length()),
                  rep("CCR7_Th",  times = CCR7_Th_specific_ligands %>% length()),
                  rep("General",  times = general_ligands %>% length())),
  ligand = c(Monocyte_specific_ligands, B_Memory_specific_ligands, B_Plasma_specific_ligands, CCR6_Th_specific_ligands, CCR7_Th_specific_ligands, general_ligands))
ligand_type_indication_df %>% head
##
active_ligand_target_links_df = nichenet_output$ligand_target_df %>% mutate(target_type = "AD_DE") %>% inner_join(ligand_type_indication_df)
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)
active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)
ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())
circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)
##
grid_col_ligand = c("General" = "lawngreen", "Monocyte_specific"="royalblue", "B_Memory_specific_specific"="darkgreen",
                    "B_Plasma"="violet","CCR6_Th_specific"="steelblue2","CCR7_Th_specific"="red")
grid_col_target = c("AD_DE" = "tomato")
grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(ligand_type = grid_col_target %>% names(), color_ligand_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand, " "))
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand, target, weight)

ligand_color = circos_links %>% distinct(ligand, color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target, color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)
grid_col = c(grid_ligand_color, grid_target_color)










