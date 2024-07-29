# NicheNets: https://github.com/saeyslab/nichenetr
install.packages("devtools")
devtools::install_github("saeyslab/nichenetr")

## Dataのdownload
## https://zenodo.org/record/3260758#.Y6jpTuzP0Ss

# Preparation:
location="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/"
ligand_target_matrix = readRDS(paste0(location,"ligand_target_matrix.rds"))
lr_network = readRDS(paste0(location,"lr_network.rds"))
weighted_networks= readRDS(paste0(location,"weighted_networks.rds"))


## Standardなmodel : muliple sender cells -> Reciever Cells
## Perform the NicheNet analysis
# indicated cell types should be cell class identities
# check via: 
# seuratObj %>% Idents() %>% table()
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seuratObj, 
  receiver = "CD8 T", 
  condition_colname = "aggregate", condition_oi = "LCMV", condition_reference = "SS", 
  sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, weighted_networks = weighted_networks, organism = "mouse")
## [1] "Read in and process NicheNet's networks"
## [1] "Define expressed ligands and receptors in receiver and sender cells"
## [1] "Perform DE analysis in receiver cell"
## [1] "Perform NicheNet ligand activity analysis"
## [1] "Infer active target genes of the prioritized ligands"
## [1] "Infer receptors of the prioritized ligands"
## [1] "Perform DE analysis in sender cells"



## One Sender to One Receiver
nichenet_output = 
  nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj,
    receiver = "CD8 T", 
    condition_colname = "aggregate", 
    condition_oi = "LCMV", 
    condition_reference = "SS", 
    sender = "DC",            # <- ここを1つにする  
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "mouse")



## One Sender to One Receiver : Autocrine
nichenet_output = 
  nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj,
    receiver = "CD8 T",          # <- receiver senderを同じ細胞にする
    condition_colname = "aggregate", 
    condition_oi = "LCMV", 
    condition_reference = "SS", 
    sender = "CD8 T",            # <- receiver senderを同じ細胞にする
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "mouse")





## All sender  to One Receiver
nichenet_output = 
  nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj,
    receiver = "CD8 T", 
    condition_colname = "aggregate", 
    condition_oi = "LCMV", 
    condition_reference = "SS", 
    sender = "all",            # <- allにすると全ての細胞との相互作用になる
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "mouse")


## Unknown sender to One Receiver
nichenet_output = 
  nichenet_seuratobj_aggregate(
    seurat_obj = seuratObj,
    receiver = "CD8 T", 
    condition_colname = "aggregate", 
    condition_oi = "LCMV", 
    condition_reference = "SS", 
    sender = "undefined",            # <- senderが不明でもreceptorのactivityからligandを推測
    ligand_target_matrix = ligand_target_matrix, 
    lr_network = lr_network, 
    weighted_networks = weighted_networks, 
    organism = "mouse")



##　応用編
## Several Receivers
receiver_celltypes_oi = c("CD4 T", "CD8 T")
# receiver_celltypes_oi = seuratObj %>% Idents() %>% unique() # for all celltypes in the dataset: use only when this would make sense biologically
nichenet_output = 
  receiver_celltypes_oi %>% 
    lapply(nichenet_seuratobj_aggregate, 
           seurat_obj = seuratObj, 
           condition_colname = "aggregate", 
           condition_oi = "LCMV", 
           condition_reference = "SS", 
           sender = c("CD4 T","Treg", "Mono", "NK", "B", "DC"), 
           ligand_target_matrix = ligand_target_matrix, 
           lr_network = lr_network, 
           weighted_networks = weighted_networks, 
           organism = "mouse")
names(nichenet_output) = receiver_celltypes_oi



## Check which ligands were top-ranked for both CD8T and CD4T and which ligands were more cell-type specific

common_ligands = intersect(nichenet_output$`CD4 T`$top_ligands, nichenet_output$`CD8 T`$top_ligands)
print("common ligands are: ")
## [1] "common ligands are: "
print(common_ligands)
##  [1] "Ebi3"   "Il15"   "Crlf2"  "H2-M3"  "App"    "Ptprc"  "Icam1"  "Ccl5"   "Cxcl10" "Tgfb1"  "Cxcl11" "Sema4d" "Cxcl9"  "H2-T23" "Cxcl16" "C3"     "Itgb1"

cd4_ligands = nichenet_output$`CD4 T`$top_ligands %>% setdiff(nichenet_output$`CD8 T`$top_ligands)
cd8_ligands = nichenet_output$`CD8 T`$top_ligands %>% setdiff(nichenet_output$`CD4 T`$top_ligands)

print("Ligands specifically regulating DE in CD4T: ")
## [1] "Ligands specifically regulating DE in CD4T: "
print(cd4_ligands)
## [1] "Cd274" "Hmgb1" "Cd28"

print("Ligands specifically regulating DE in CD8T: ")
## [1] "Ligands specifically regulating DE in CD8T: "
print(cd8_ligands)
## [1] "Adam17" "Anxa1"  "Sell"



### Condition x Cell typeなどで，同一細胞の状態が異なる間での相互作用などを解析
## Seuratで新たなmeta dataを作成
seuratObj@meta.data$celltype = 
  paste(seuratObj@meta.data$celltype,
        seuratObj@meta.data$aggregate, 
        sep = "_")

  seuratObj@meta.data$celltype %>% table()
## .
##     B_LCMV       B_SS CD4 T_LCMV   CD4 T_SS CD8 T_LCMV   CD8 T_SS    DC_LCMV      DC_SS  Mono_LCMV    Mono_SS    NK_LCMV      NK_SS  Treg_LCMV    Treg_SS 
##        344         38       1961        601       1252        393         14          4         75         15         94         37        146         53

  seuratObj = SetIdent(seuratObj,value = "celltype")
  
## 作成した新たなmeta dataで実施
nichenet_output = nichenet_seuratobj_cluster_de(
  seurat_obj = seuratObj, 
  receiver_reference = "CD8 T_SS", 
  receiver_affected = "CD8 T_LCMV", 
  sender = c("DC_LCMV","Mono_LCMV"), 
  ligand_target_matrix = ligand_target_matrix, 
  lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  organism = "mouse")










############ 結果の確認 ################ 
#To get a list of the 20 top-ranked ligands: run the following command
nichenet_output$top_ligands
nichenet_output$ligand_expression_dotplot
nichenet_output$ligand_differential_expression_heatmap
nichenet_output$ligand_target_heatmap
nichenet_output$ligand_target_heatmap + 
  scale_fill_gradient2(low = "whitesmoke",  high = "royalblue", breaks = c(0,0.0045,0.009)) + 
  xlab("Recivor response genes") + 
  ylab("Prioritized sender genes")
nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
nichenet_output$ligand_target_df # weight column = regulatory potential
# To get a list of the top-predicted target genes of the 20 top-ranked ligands
nichenet_output$top_targets

## Seurat objectを使って Targetの発現を確認
## Top20 ligands for dot plot 
DotPlot(wh.int %>% subset(idents = "Endothelial"), 
        features = nichenet_output$top_targets %>% rev(), 
        split.by = "State") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(colour = "black", size= 14, face="bold"),
        axis.text.y = element_text(colour = "black", size= 14, face="bold"))

## Vlnplot
VlnPlot(wh.int %>% subset(idents = "Endothelial"), 
        features = c("Robo4"), 
        split.by = "State", pt.size = 0, combine = FALSE)


## [[1]]
# To visualize ligand activities, expression, differential expression and target genes of ligands
## one of the most important summary figures of the NicheNet analysis.
nichenet_output$ligand_activity_target_heatmap


## Inferred ligand-receptor interactions for top-ranked ligands
nichenet_output$ligand_receptor_heatmap
## data file 
  nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
  nichenet_output$ligand_receptor_df # weight column accords to number of data sources that document this interaction
  
# To get a list of the receptors of the 20 top-ranked ligands
nichenet_output$top_receptors
# Top20 ligands for dot plot 
DotPlot(wh.int %>% subset(idents = "Endothelial"),
        features = nichenet_output$top_receptors %>% rev(), 
        split.by = "State") + 
  RotatedAxis() +
  theme(axis.text.x = element_text(colour = "black", size= 14, face="bold"),
        axis.text.y = element_text(colour = "black", size= 14, face="bold"))


## NicheNet also infers the receiver cell receptors of these top-ranked ligands. You can run following command for a heatmap visualization of the ligand-receptor links:
nichenet_output$ligand_receptor_heatmap_bonafide
nichenet_output$ligand_receptor_matrix %>% .[1:10,1:6]
nichenet_output$ligand_receptor_matrix_bonafide
nichenet_output$ligand_receptor_df_bonafide


nichenet_output$geneset_oi
nichenet_output$background_expressed_genes %>% length()






