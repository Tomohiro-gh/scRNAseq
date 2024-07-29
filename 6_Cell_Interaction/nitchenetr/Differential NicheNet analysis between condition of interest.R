## Differential NicheNet analysis between conditions of interest
## https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet_pEMT.md

## ２つのsampleの状態での相互作用を比べる

library(nichenetr) #v1.1.1 at 01/21/23
library(RColorBrewer)
library(tidyverse)
library(Seurat) 

## 0 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Data preparation

seurat_obj = readRDS("/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/hnscc_expression.rds")

DimPlot(seurat_obj, group.by = "celltype") # user adaptation required on own dataset
DimPlot(seurat_obj, group.by = "pEMT") # user adaptation required on own dataset

table(seurat_obj@meta.data$celltype, seurat_obj@meta.data$pEMT) # cell types vs conditions # user adaptation required on own dataset

seurat_obj@meta.data$celltype_aggregate = 
  paste(
    seurat_obj@meta.data$celltype, 
    seurat_obj@meta.data$pEMT,
    sep = "_") # user adaptation required on own dataset
DimPlot(seurat_obj, group.by = "celltype_aggregate")

seurat_obj@meta.data$celltype_aggregate %>% 
  table() %>% 
  sort(decreasing = TRUE)

celltype_id = "celltype_aggregate" # metadata column name of the cell type of interest
seurat_obj = 
  SetIdent(seurat_obj, 
           value = seurat_obj[[celltype_id]])

### Read in the NicheNet ligand-receptor network and ligand-target matrix
location="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/"
ligand_target_matrix = 
  readRDS(paste0(location,"ligand_target_matrix.rds"))

lr_network = readRDS(paste0(location,"lr_network.rds"))
lr_network = 
  lr_network %>% 
  mutate(bonafide = ! database %in% 
           c("ppi_prediction","ppi_prediction_go"))
lr_network = 
  lr_network %>% 
  dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor, bonafide)

organism = "human" # user adaptation required on own dataset

if(organism == "mouse"){
  lr_network = 
    lr_network %>% 
    mutate(ligand = convert_human_to_mouse_symbols(ligand), 
           receptor = convert_human_to_mouse_symbols(receptor)) %>%
    drop_na()
  
  colnames(ligand_target_matrix) = 
    ligand_target_matrix %>% 
    colnames() %>% 
    convert_human_to_mouse_symbols()
  
  rownames(ligand_target_matrix) = 
    ligand_target_matrix %>% 
    rownames() %>% 
    convert_human_to_mouse_symbols()
  ligand_target_matrix = 
    ligand_target_matrix %>% 
    .[!is.na(rownames(ligand_target_matrix)), 
      !is.na(colnames(ligand_target_matrix))]
}

## 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Define the niches/microenvironments of interest
niches = list(
  "pEMT_High_niche" = list(
    "sender" = c("myofibroblast_High", "Endothelial_High", "CAF_High", "T.cell_High", "Myeloid_High"),
    "receiver" = c("Malignant_High")),
  "pEMT_Low_niche" = list(
    "sender" = c("myofibroblast_Low",  "Endothelial_Low", "CAF_Low"),
    "receiver" = c("Malignant_Low"))
) # user adaptation required on own dataset



##  2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate differential expression between the niches
assay_oi = "SCT" # other possibilities: RNA,...

DE_sender = 
  calculate_niche_de(
    seurat_obj = 
      seurat_obj %>% 
      subset(features = lr_network$ligand %>% 
               unique()), 
    niches = niches, 
    type = "sender", 
    assay_oi = assay_oi) # only ligands important for sender cell types

DE_receiver = 
  calculate_niche_de(
    seurat_obj = 
      seurat_obj %>% 
      subset(features = lr_network$receptor %>% 
               unique()), 
    niches = niches, 
    type = "receiver", 
    assay_oi = assay_oi) # only receptors now, later on: DE analysis to find targets

DE_sender = 
  DE_sender %>% 
  mutate(
    avg_log2FC = 
      ifelse(avg_log2FC == Inf, 
             max(avg_log2FC[is.finite(avg_log2FC)]), 
             ifelse(avg_log2FC == -Inf, 
                    min(avg_log2FC[is.finite(avg_log2FC)]), 
                    avg_log2FC)))
DE_receiver = 
  DE_receiver %>%
  mutate(
    avg_log2FC = 
      ifelse(avg_log2FC == Inf, 
             max(avg_log2FC[is.finite(avg_log2FC)]), 
             ifelse(avg_log2FC == -Inf, 
                    min(avg_log2FC[is.finite(avg_log2FC)]), 
                    avg_log2FC)))

## Process DE results:
expression_pct = 0.10
DE_sender_processed = 
  process_niche_de(
    DE_table = DE_sender, 
    niches = niches, 
    expression_pct = expression_pct, 
    type = "sender")
DE_receiver_processed = 
  process_niche_de(
    DE_table = DE_receiver, 
    niches = niches, 
    expression_pct = expression_pct, 
    type = "receiver")

## Combine sender-receiver DE based on L-R pairs:
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = 
  combine_sender_receiver_de(
    DE_sender_processed, 
    DE_receiver_processed, 
    lr_network, 
    specificity_score = 
      specificity_score_LR_pairs)  

## 3 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Optional: Calculate differential expression between the different spatial regions

include_spatial_info_sender = FALSE # if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
## ここは　skip

spatial_info = 
  tibble(
    celltype_region_oi = "CAF_High", 
    celltype_other_region = "myofibroblast_High", 
    niche =  "pEMT_High_niche", 
    celltype_type = "sender") # user adaptation required on own dataset
specificity_score_spatial = "lfc"

# this is how this should be defined if you don't have spatial info
# mock spatial info
if(
  include_spatial_info_sender == FALSE & 
  include_spatial_info_receiver == FALSE){
  spatial_info = 
    tibble(
      celltype_region_oi = NA, 
      celltype_other_region = NA) %>% 
    mutate(niche =  niches %>% 
             names() %>% 
             head(1), 
           celltype_type = "sender")
} 

## TRUEの場合は，こっち
if(
  include_spatial_info_sender == TRUE){
  sender_spatial_DE = 
    calculate_spatial_DE(
      seurat_obj = 
        seurat_obj %>% 
        subset(features = lr_network$ligand %>% 
                 unique()), 
      spatial_info = spatial_info %>% 
        filter(celltype_type == "sender"))
  sender_spatial_DE_processed = 
    process_spatial_de(
      DE_table = sender_spatial_DE, 
      type = "sender", 
      lr_network = lr_network, 
      expression_pct = expression_pct, 
      specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for sender celltypes in which the spatial is not known / not of importance
  sender_spatial_DE_others = 
    get_non_spatial_de(
      niches = niches, 
      spatial_info = spatial_info, 
      type = "sender", 
      lr_network = lr_network)
  sender_spatial_DE_processed = 
    sender_spatial_DE_processed %>% 
    bind_rows(sender_spatial_DE_others)
  
  sender_spatial_DE_processed = 
    sender_spatial_DE_processed %>% 
    mutate(
      scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
  
} else {
  # # add a neutral spatial score for all sender celltypes (for none of them, spatial is relevant in this case)
  sender_spatial_DE_processed = 
    get_non_spatial_de(
      niches = niches, 
      spatial_info = spatial_info, 
      type = "sender", 
      lr_network = lr_network)
  sender_spatial_DE_processed = 
    sender_spatial_DE_processed %>% 
    mutate(
      scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))  
  
}
## [1] "Calculate Spatial DE between: CAF_High and myofibroblast_High"

## Recieverの方の処理
if(
  include_spatial_info_receiver == TRUE){
  receiver_spatial_DE = 
    calculate_spatial_DE(
      seurat_obj = seurat_obj %>% 
        subset(features = lr_network$receptor %>% 
                 unique()), 
      spatial_info = spatial_info %>% 
        filter(celltype_type == "receiver"))
  receiver_spatial_DE_processed = 
    process_spatial_de(
      DE_table = receiver_spatial_DE, 
      type = "receiver", 
      lr_network = lr_network, 
      expression_pct = expression_pct, 
      specificity_score = specificity_score_spatial)
  
  # add a neutral spatial score for receiver celltypes in which the spatial is not known / not of importance
  receiver_spatial_DE_others = 
    get_non_spatial_de(
      niches = niches, 
      spatial_info = spatial_info, 
      type = "receiver", 
      lr_network = lr_network)
  receiver_spatial_DE_processed = 
    receiver_spatial_DE_processed %>% 
    bind_rows(receiver_spatial_DE_others)
  
  receiver_spatial_DE_processed =
    receiver_spatial_DE_processed %>% 
    mutate(
      scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
  
} else {
  # # add a neutral spatial score for all receiver celltypes (for none of them, spatial is relevant in this case)
  receiver_spatial_DE_processed = 
    get_non_spatial_de(
      niches = niches, 
      spatial_info = spatial_info, 
      type = "receiver", 
      lr_network = lr_network)
  receiver_spatial_DE_processed = 
    receiver_spatial_DE_processed %>% 
    mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}



## 4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate ligand activities and infer active ligand-target links
lfc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

DE_receiver_targets = 
  calculate_niche_de_targets(
    seurat_obj = seurat_obj, 
    niches = niches, 
    lfc_cutoff = lfc_cutoff, 
    expression_pct = expression_pct, 
    assay_oi = assay_oi) 
## [1] "Calculate receiver DE between: Malignant_High and Malignant_Low"
## [1] "Calculate receiver DE between: Malignant_Low and Malignant_High"
DE_receiver_processed_targets = 
  process_receiver_target_de(
    DE_receiver_targets = DE_receiver_targets, 
    niches = niches, 
    expression_pct = expression_pct, 
    specificity_score = specificity_score_targets)

background = 
  DE_receiver_processed_targets  %>% 
  pull(target) %>% 
  unique()

geneset_niche1 = 
  DE_receiver_processed_targets %>% 
  filter(
    receiver == niches[[1]]$receiver & 
      target_score >= lfc_cutoff & 
      target_significant == 1 & 
      target_present == 1) %>% 
  pull(target) %>% 
  unique()

geneset_niche2 = 
  DE_receiver_processed_targets %>% 
  filter(
    receiver == niches[[2]]$receiver & 
      target_score >= lfc_cutoff & 
      target_significant == 1 & 
    arget_present == 1) %>% 
  pull(target) %>% 
  unique()

geneset_niche1 %>% 
  setdiff(rownames(ligand_target_matrix))

geneset_niche2 %>% 
  setdiff(rownames(ligand_target_matrix))


length(geneset_niche1)
## [1] 1668
length(geneset_niche2)
## [1] 2889


## recommend having between 20 and 1000 genes in the geneset of interest, and a background of at least 5000 genes for a proper ligand activity analysis. 
##  We recommend using a cutoff of 0.15 if you have > 2 receiver cells/niches to compare and use the min_lfc as specificity score. If you have only 2 receivers/niche, we recommend using a higher threshold (such as using 0.25). 

#### 変動遺伝子が多く出過ぎたら，lfc_cutoff = 0.75 の値を高くしてみる
lfc_cutoff = 0.75 
specificity_score_targets = "min_lfc"

DE_receiver_processed_targets = 
  process_receiver_target_de(
    DE_receiver_targets = DE_receiver_targets, 
    niches = niches, 
    expression_pct = expression_pct, 
    specificity_score = specificity_score_targets)

background = 
  DE_receiver_processed_targets  %>% 
  pull(target) %>% 
  unique()

geneset_niche1 = 
  DE_receiver_processed_targets %>% 
  filter(
    receiver == niches[[1]]$receiver & 
      target_score >= lfc_cutoff & 
      target_significant == 1 & 
      target_present == 1) %>% 
  pull(target) %>% 
  unique()

geneset_niche2 = 
  DE_receiver_processed_targets %>% 
  filter(
    receiver == niches[[2]]$receiver & 
      target_score >= lfc_cutoff & 
      target_significant == 1 & 
      target_present == 1) %>% 
  pull(target) %>% 
  unique()
###  optional 終わり<<<<<<<<<<<<<<<< 


top_n_target = 250

niche_geneset_list = list(
  "pEMT_High_niche" = list(
    "receiver" = niches[[1]]$receiver,
    "geneset" = geneset_niche1,
    "background" = background),
  "pEMT_Low_niche" = list(
    "receiver" = niches[[2]]$receiver,
    "geneset" = geneset_niche2 ,
    "background" = background)
)


ligand_activities_targets = 
  get_ligand_activities_targets(
    niche_geneset_list = niche_geneset_list, 
    ligand_target_matrix = ligand_target_matrix, 
    top_n_target = top_n_target)


## 5 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)

# In this step, we will calculate average (scaled) expression, and fraction of expression, of ligands, receptors, and target genes across all cell types of interest. Now this is here demonstrated via the DotPlot function of Seurat, but this can also be done via other ways of course.

features_oi = 
  union(lr_network$ligand, lr_network$receptor) %>% 
  union(ligand_activities_targets$target) %>% 
  setdiff(NA)

dotplot = 
  suppressWarnings(
    Seurat::DotPlot(seurat_obj %>% 
                      subset(idents = niches %>% 
                               unlist() %>% 
                               unique()), 
                    features = features_oi, 
                    assay = assay_oi))
exprs_tbl = 
  dotplot$data %>% 
  as_tibble()

exprs_tbl = 
  exprs_tbl %>% 
  rename(celltype = id, 
         gene = features.plot, 
         expression = avg.exp, 
         expression_scaled = avg.exp.scaled, 
         fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% 
  as_tibble() %>% 
  select(celltype, gene, expression, expression_scaled, fraction) %>% 
  distinct() %>% 
  arrange(gene) %>% 
  mutate(gene = as.character(gene))

exprs_tbl_ligand = 
  exprs_tbl %>% 
  filter(gene %in% lr_network$ligand) %>% 
  rename(sender = celltype, 
         ligand = gene, 
         ligand_expression = expression, 
         ligand_expression_scaled = expression_scaled, 
         ligand_fraction = fraction) 

exprs_tbl_receptor = 
  exprs_tbl %>% 
  filter(gene %in% lr_network$receptor) %>% 
  rename(receiver = celltype, 
         receptor = gene, 
         receptor_expression = expression, 
         receptor_expression_scaled = expression_scaled,
         receptor_fraction = fraction)

exprs_tbl_target = 
  exprs_tbl %>% 
  filter(gene %in% ligand_activities_targets$target) %>% 
  rename(receiver = celltype, 
         target = gene, 
         target_expression = expression, 
         target_expression_scaled = expression_scaled, 
         target_fraction = fraction)

## 6 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Expression fraction and receptor 
exprs_sender_receiver = lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>% 
  inner_join(DE_sender_receiver %>% 
               distinct(niche, sender, receiver))

## >>>>> Different points <<<<<<
ligand_scaled_receptor_expression_fraction_df = 
  exprs_sender_receiver %>%
  group_by(ligand, receiver) %>% 
  mutate(
    rank_receptor_expression = dense_rank(receptor_expression), 
    rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% 
  mutate(
    ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% 
  distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>% 
  distinct() %>% 
  ungroup() 

## 7 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Prioritization of ligand-receptor and ligand-target links

##In this step, we will combine all the above calculated information to prioritize ligand-receptor-target links. We scale every property of interest between 0 and 1, and the final prioritization score is a weighted sum of the scaled scores of all the properties of interest.
## 優先順位をつける

prioritizing_weights = 
  c("scaled_ligand_score" = 5,
    "scaled_ligand_expression_scaled" = 1,
    "ligand_fraction" = 1,
    "scaled_ligand_score_spatial" = 2,
    "scaled_receptor_score" = 0.5,
    "scaled_receptor_expression_scaled" = 0.5,
    "receptor_fraction" = 1,
    "ligand_scaled_receptor_expression_fraction" = 1,
    "scaled_receptor_score_spatial" = 0,
    "scaled_activity" = 0,
    "scaled_activity_normalized" = 1,
    "bona_fide" = 1)

## Outputはリストにして
output = list(
  DE_sender_receiver = DE_sender_receiver, 
  ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
  #sender_spatial_DE_processed = sender_spatial_DE_processed, 
  #receiver_spatial_DE_processed = receiver_spatial_DE_processed,
  ligand_activities_targets = ligand_activities_targets, 
  DE_receiver_processed_targets = DE_receiver_processed_targets, 
  exprs_tbl_ligand = exprs_tbl_ligand,  
  exprs_tbl_receptor = exprs_tbl_receptor, 
  exprs_tbl_target = exprs_tbl_target)

prioritization_tables = 
  get_prioritization_tables(output, prioritizing_weights)



prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[1]]$receiver) %>% 
  head(10)

prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[1]]$receiver) %>% 
  head(10)

prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[2]]$receiver) %>% 
  head(10)

prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[2]]$receiver) %>%
  head(10)

prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(receiver == niches[[3]]$receiver) %>% 
  head(10)

prioritization_tables$prioritization_tbl_ligand_target %>% 
  filter(receiver == niches[[3]]$receiver) %>% 
  head(10)

prioritization_tables$prioritization_tbl_ligand_receptor = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  mutate(receiver = factor(receiver, 
                           levels = c("KCs","MoMac1","MoMac2")), 
         niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac2_niche"))) 

prioritization_tables$prioritization_tbl_ligand_target = 
  prioritization_tables$prioritization_tbl_ligand_target %>% 
  mutate(receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), 
         niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac2_niche"))) 


## 8 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Visualization of the Differential NicheNet output

# Before visualization, we need to define the most important ligand-receptor pairs per niche. 
# We will do this by first determining for which niche the highest score is found for each ligand/ligand-receptor pair. And then getting the top 50 ligands per niche.
top_ligand_niche_df = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, 
         sender, 
         receiver, 
         ligand, 
         receptor, 
         prioritization_score) %>% 
  group_by(ligand) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand, receptor, niche) %>% 
  rename(top_niche = niche)

top_ligand_receptor_niche_df =
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, 
         sender, 
         receiver, 
         ligand, 
         receptor, 
         prioritization_score) %>% 
  group_by(ligand, receptor) %>%
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand, receptor, niche) %>%
  rename(top_niche = niche)

ligand_prioritized_tbl_oi =
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, 
         sender, 
         receiver, 
         ligand, 
         prioritization_score) %>%
  group_by(ligand, niche) %>% 
  top_n(1, prioritization_score) %>%
  ungroup() %>%
  distinct() %>%
  inner_join(top_ligand_niche_df) %>%
  filter(niche == top_niche) %>%
  group_by(niche) %>%
  top_n(50, prioritization_score) %>%
  ungroup() 

# Now we will look first at the top ligand-receptor pairs for KCs (here, we will take the top 2 scoring receptors per prioritized ligand)
receiver_oi = "Malignant_High" 

filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  pull(ligand) %>% 
  unique()

prioritized_tbl_oi = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, 
         sender, 
         receiver, 
         ligand, 
         receptor, 
         ligand_receptor, 
         prioritization_score) %>% 
  distinct() %>% 
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% 
  ungroup() 


## Visualization: minimum LFC compared to other niches
lfc_plot = 
  make_ligand_receptor_lfc_plot(
    receiver_oi, 
    prioritized_tbl_oi, 
    prioritization_tables$prioritization_tbl_ligand_receptor,
    plot_legend = FALSE,
    heights = NULL, 
    widths = NULL)
lfc_plot

#Show the spatialDE as additional information
lfc_plot_spatial = 
  make_ligand_receptor_lfc_spatial_plot(
    receiver_oi, 
    prioritized_tbl_oi, 
    prioritization_tables$prioritization_tbl_ligand_receptor, 
    ligand_spatial = include_spatial_info_sender, 
    receptor_spatial = include_spatial_info_receiver, 
    plot_legend = FALSE, 
    heights = NULL, 
    widths = NULL)
lfc_plot_spatial



## Ligand expression, activity and target genes
## Active target gene inference - cf Default NicheNet
exprs_activity_target_plot =
  make_ligand_activity_target_exprs_plot(
    receiver_oi, 
    prioritized_tbl_oi,  
    prioritization_tables$prioritization_tbl_ligand_receptor, 
    prioritization_tables$prioritization_tbl_ligand_target,
    output$exprs_tbl_ligand, 
    output$exprs_tbl_target,
    lfc_cutoff,
    ligand_target_matrix,
    plot_legend = FALSE, 
    heights = NULL,
    widths = NULL)
exprs_activity_target_plot$combined_plot
## important: ligand-receptor pairs with both high differential expression (or condition-specificity) and ligand activity (=target gene enrichment) are very interesting predictions as key regulators of your intercellular communication process of interest !
  
  
#If this plot contains too much information because we look at many hits (top 50 ligands), you can make this plot of course for less ligands as well, eg for the top20.
filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>%
  top_n(20, prioritization_score) %>% 
  pull(ligand) %>%
  unique()

prioritized_tbl_oi = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(ligand %in% filtered_ligands) %>% 
  select(niche,
         sender,
         receiver,
         ligand, 
         receptor, 
         ligand_receptor,
         prioritization_score) %>% 
  distinct() %>% 
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>%
  filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% 
  ungroup() 

exprs_activity_target_plot = 
  make_ligand_activity_target_exprs_plot(
    receiver_oi, 
    prioritized_tbl_oi,  
    prioritization_tables$prioritization_tbl_ligand_receptor, 
    prioritization_tables$prioritization_tbl_ligand_target,
    output$exprs_tbl_ligand, 
    output$exprs_tbl_target,
    lfc_cutoff,
    ligand_target_matrix,
    plot_legend = FALSE, 
    heights = NULL,
    widths = NULL)
exprs_activity_target_plot$combined_plot




## Circos plot of prioritized ligand-receptor pairs
# Because a top50 is too much to visualize in a circos plot, we will only visualize the top 15.
filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(15, prioritization_score) %>% 
  pull(ligand) %>% 
  unique()

prioritized_tbl_oi =
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, 
         sender,
         receiver,
         ligand, 
         receptor,
         ligand_receptor,
         prioritization_score) %>% 
  distinct() %>% 
  inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% 
  ungroup() 

colors_sender = 
  brewer.pal(
    n = prioritized_tbl_oi$sender %>% 
      unique() %>% 
      sort() %>% 
      length(),
    name = 'Spectral') %>% 
  magrittr::set_names(
    prioritized_tbl_oi$sender %>% 
      unique() %>% 
      sort())
colors_receiver = 
  c("lavender")  %>% 
  magrittr::set_names(
    prioritized_tbl_oi$receiver %>% 
      unique() %>% 
      sort())

circos_output = 
  make_circos_lr(
    prioritized_tbl_oi, 
    colors_sender, 
    colors_receiver)

## Visualization for the other condition: pEMT-low
receiver_oi = "Malignant_Low"  
filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(50, prioritization_score) %>%
  pull(ligand) %>% 
  unique()

prioritized_tbl_oi = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, 
         sender,
         receiver,
         ligand, 
         receptor,
         ligand_receptor, 
         prioritization_score) %>%
  distinct() %>%
  inner_join(top_ligand_receptor_niche_df) %>%
  group_by(ligand) %>%
  filter(receiver == receiver_oi) %>%
  top_n(2, prioritization_score) %>% 
  ungroup() 

lfc_plot = 
  make_ligand_receptor_lfc_plot(
    receiver_oi,
    prioritized_tbl_oi, 
    prioritization_tables$prioritization_tbl_ligand_receptor,
    plot_legend = FALSE,
    heights = NULL, 
    widths = NULL)
lfc_plot