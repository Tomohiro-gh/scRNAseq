## Differential NicheNet analysis between conditions of interest
## https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet.md


library(nichenetr) #v1.1.1 at 01/21/23
library(RColorBrewer)
library(tidyverse)
library(Seurat) 


# deposit data:
# https://zenodo.org/record/5840787/#.Y8t6BOzP0eU


seurat_obj =
  readRDS("/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/seurat_obj_subset_integrated_zonation.rds")

#We will now also check the number of cells per cell type condition combination
table(seurat_obj@meta.data$celltype, seurat_obj@meta.data$pEMT)


seurat_obj = SetIdent(seurat_obj, value = "celltype")


## 0 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Read in the NicheNet ligand-receptor network and ligand-target matrix
location="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/"
ligand_target_matrix = readRDS(paste0(location,"ligand_target_matrix.rds"))

lr_network = readRDS(paste0(location,"lr_network.rds"))
lr_network = 
  lr_network %>% 
  mutate(bonafide = ! database %in% 
           c("ppi_prediction","ppi_prediction_go"))
lr_network = 
  lr_network %>% 
  dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor, bonafide)

weighted_networks = readRDS(paste0(location,"weighted_networks.rds"))

## Mouseの場合，humanへ変換する
organism = "mouse" 

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

## ! Important: your receiver cell type should consist of 1 cluster!

niches = list(
  "KC_niche" = list(
    "sender" = c("LSECs_portal","Hepatocytes_portal","Stellate cells_portal"),
    "receiver" = c("KCs")),
  "MoMac2_niche" = list(
    "sender" = c("Cholangiocytes","Fibroblast 2"),
    "receiver" = c("MoMac2")),
  "MoMac1_niche" = list(
    "sender" = c("Capsule fibroblasts","Mesothelial cells"),
    "receiver" = c("MoMac1"))
)

## 2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate differential expression between the niches
## Calculate DE
assay_oi = "SCT" # other possibilities: RNA,...
seurat_obj = 
  PrepSCTFindMarkers(
    seurat_obj, 
    assay = "SCT", 
    verbose = FALSE) # RNAならこれは不要

DE_sender = 
  calculate_niche_de(
    seurat_obj = 
      seurat_obj %>% 
      subset(
        features = lr_network$ligand %>% 
          intersect(rownames(seurat_obj))), 
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
  mutate(avg_log2FC = 
           ifelse(avg_log2FC == Inf,
                  max(avg_log2FC[is.finite(avg_log2FC)]),
                  ifelse(avg_log2FC == -Inf,
           min(avg_log2FC[is.finite(avg_log2FC)]), 
           avg_log2FC)))

#Process DE results:
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

# Combine sender-receiver DE based on L-R pairs:
specificity_score_LR_pairs = "min_lfc"
DE_sender_receiver = 
  combine_sender_receiver_de(
    DE_sender_processed, 
    DE_receiver_processed, 
    lr_network, 
    specificity_score = specificity_score_LR_pairs)



## 3 >>>>>>>>>>>>>>>>>>>>>>>>>>d>>>>>>>>>>>>
## Optional: Calculate differential expression between the different spatial regions

## ここは省略

## 4 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate ligand activities and infer active ligand-target links
fc_cutoff = 0.15 # recommended for 10x as min_lfc cutoff. 
specificity_score_targets = "min_lfc"

DE_receiver_targets = 
  calculate_niche_de_targets(
    seurat_obj = seurat_obj, 
    niches = niches, 
    lfc_cutoff = lfc_cutoff, 
    expression_pct = expression_pct, 
    assay_oi = assay_oi) 
## [1] "Calculate receiver DE between: KCs and MoMac2" "Calculate receiver DE between: KCs and MoMac1"
## [1] "Calculate receiver DE between: MoMac2 and KCs"    "Calculate receiver DE between: MoMac2 and MoMac1"
## [1] "Calculate receiver DE between: MoMac1 and KCs"    "Calculate receiver DE between: MoMac1 and MoMac2"
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
geneset_KC = 
  DE_receiver_processed_targets %>% 
  filter(receiver == niches$KC_niche$receiver & 
           target_score >= lfc_cutoff & 
           target_significant == 1 & 
           target_present == 1) %>%
  pull(target) %>% 
  unique()

geneset_MoMac2 = 
  DE_receiver_processed_targets %>% 
  filter(receiver == niches$MoMac2_niche$receiver & 
           target_score >= lfc_cutoff & 
           target_significant == 1 & 
           target_present == 1) %>% 
  pull(target) %>% 
  unique()

geneset_MoMac1 = 
  DE_receiver_processed_targets %>% 
  filter(receiver == niches$MoMac1_niche$receiver & 
           target_score >= lfc_cutoff & 
           target_significant == 1 & 
           target_present == 1) %>% 
  pull(target) %>% 
  unique()



top_n_target = 250
## 
niche_geneset_list = list(
  "KC_niche" = list(
    "receiver" = "KCs",
    "geneset" = geneset_KC,
    "background" = background),
  "MoMac1_niche" = list(
    "receiver" = "MoMac1",
    "geneset" = geneset_MoMac1 ,
    "background" = background),
  "MoMac2_niche" = list(
    "receiver" = "MoMac2",
    "geneset" = geneset_MoMac2 ,
    "background" = background)  
)

ligand_activities_targets = 
  get_ligand_activities_targets(
    niche_geneset_list = niche_geneset_list, 
    ligand_target_matrix = ligand_target_matrix, 
    top_n_target = top_n_target)




## 5 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate (scaled) expression of ligands, receptors and targets across cell types of interest (log expression values and expression fractions)
features_oi = 
  union(lr_network$ligand, 
        lr_network$receptor) %>% 
  union(ligand_activities_targets$target) %>% 
  setdiff(NA)

dotplot = 
  suppressWarnings(
    Seurat::DotPlot(
      seurat_obj %>% 
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
  rename(
    celltype = id, 
    gene = features.plot, 
    expression = avg.exp, 
    expression_scaled = avg.exp.scaled, 
    fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% 
  as_tibble() %>% 
  select(
    celltype, 
    gene, 
    expression, 
    expression_scaled, 
    fraction) %>% 
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

exprs_tbl_ligand = 
  exprs_tbl_ligand %>%  
  mutate(
    scaled_ligand_expression_scaled = 
      scale_quantile_adapted(ligand_expression_scaled)) %>%
  mutate(ligand_fraction_adapted = ligand_fraction) %>%
  mutate_cond(ligand_fraction >= expression_pct,
              ligand_fraction_adapted = expression_pct) %>%
  mutate(scaled_ligand_fraction_adapted =
           scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = 
  exprs_tbl_receptor %>% 
  mutate(scaled_receptor_expression_scaled =
           scale_quantile_adapted(receptor_expression_scaled)) %>% 
  mutate(receptor_fraction_adapted = 
           receptor_fraction) %>% 
  mutate_cond(receptor_fraction >= expression_pct,
              receptor_fraction_adapted = expression_pct) %>% 
  mutate(scaled_receptor_fraction_adapted =
           scale_quantile_adapted(receptor_fraction_adapted))

## 6 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Expression fraction and receptor
exprs_sender_receiver = 
  lr_network %>% 
  inner_join(exprs_tbl_ligand, by = c("ligand")) %>%
  inner_join(exprs_tbl_receptor, by = c("receptor")) %>%
  inner_join(DE_sender_receiver %>% 
               distinct(niche, sender, receiver))

ligand_scaled_receptor_expression_fraction_df =
  exprs_sender_receiver %>% 
  group_by(ligand, receiver) %>% 
  mutate(rank_receptor_expression = dense_rank(receptor_expression),
         rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% 
  mutate(ligand_scaled_receptor_expression_fraction = 
           0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% 
  distinct(ligand, 
           receptor, 
           receiver,
           ligand_scaled_receptor_expression_fraction,
           bonafide) %>% 
  distinct() %>% 
  ungroup() 



## 7 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Prioritization of ligand-receptor and ligand-target links

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

output = 
  list(DE_sender_receiver = DE_sender_receiver,
       ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
       sender_spatial_DE_processed = sender_spatial_DE_processed, 
       receiver_spatial_DE_processed = receiver_spatial_DE_processed,
       ligand_activities_targets = ligand_activities_targets,
       DE_receiver_processed_targets = DE_receiver_processed_targets,
       exprs_tbl_ligand = exprs_tbl_ligand,
       exprs_tbl_receptor = exprs_tbl_receptor,
       exprs_tbl_target = exprs_tbl_target)

prioritization_tables = 
  get_prioritization_tables(
    output, prioritizing_weights)

# 確認
prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(receiver == niches[[1]]$receiver) %>% 
  head(10)
# 確認
prioritization_tables$prioritization_tbl_ligand_target %>%
  filter(receiver == niches[[1]]$receiver) %>% 
  head(10)
# 確認
prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(receiver == niches[[2]]$receiver) %>% 
  head(10)
# 確認
prioritization_tables$prioritization_tbl_ligand_target %>%
  filter(receiver == niches[[2]]$receiver) %>% 
  head(10)
# 確認
prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(receiver == niches[[3]]$receiver) %>% 
  head(10)
# 確認
prioritization_tables$prioritization_tbl_ligand_target %>%
  filter(receiver == niches[[3]]$receiver) %>% 
  head(10)


prioritization_tables$prioritization_tbl_ligand_receptor =
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  mutate(
    receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), 
    niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac2_niche"))) 
prioritization_tables$prioritization_tbl_ligand_target =
  prioritization_tables$prioritization_tbl_ligand_target %>%
  mutate(
    receiver = factor(receiver, levels = c("KCs","MoMac1","MoMac2")), 
    niche = factor(niche, levels = c("KC_niche","MoMac1_niche","MoMac2_niche"))) 



##  8 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Visualization of the Differential NicheNet output
top_ligand_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand, receptor, niche) %>% 
  rename(top_niche = niche)

top_ligand_receptor_niche_df = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand, receptor) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand, receptor, niche) %>% 
  rename(top_niche = niche)

ligand_prioritized_tbl_oi =
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, prioritization_score) %>% 
  group_by(ligand, niche) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  distinct() %>% 
  inner_join(top_ligand_niche_df) %>% 
  filter(niche == top_niche) %>% 
  group_by(niche) %>% 
  top_n(50, prioritization_score) %>% 
  ungroup() # get the top50 ligands per niche


# Now we will look first at the top ligand-receptor pairs for KCs (here, we will take the top 2 scoring receptors per prioritized ligand)
receiver_oi = "KCs" 

filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  pull(ligand) %>% 
  unique()

prioritized_tbl_oi = prioritization_tables$prioritization_tbl_ligand_receptor %>%
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% 
  inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% 
  ungroup() 

# Visualization: minimum LFC compared to other niches
lfc_plot = 
  make_ligand_receptor_lfc_plot(
    receiver_oi, prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor, 
    plot_legend = FALSE, 
    heights = NULL, 
    widths = NULL)
lfc_plot

## Show the spatialDE as additional information
## optionalで指定した時のみ
lfc_plot_spatial = 
  make_ligand_receptor_lfc_spatial_plot(
    receiver_oi, prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor, 
    ligand_spatial = include_spatial_info_sender,
    receptor_spatial = include_spatial_info_receiver,
    plot_legend = FALSE,
    heights = NULL,
    widths = NULL)
lfc_plot_spatial


########
###Ligand expression, activity and target genes
exprs_activity_target_plot =
  make_ligand_activity_target_exprs_plot(
    receiver_oi, prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor,
    prioritization_tables$prioritization_tbl_ligand_target,
    output$exprs_tbl_ligand,
    output$exprs_tbl_target, lfc_cutoff,
    ligand_target_matrix,
    plot_legend = FALSE,
    heights = NULL,
    widths = NULL)
exprs_activity_target_plot$combined_plot

#If this plot contains too much information because we look at many hits (top 50 ligands), you can make this plot of course for less ligands as well, eg for the top20.
# top 20に減らす場合
filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(15, prioritization_score) %>% 
  pull(ligand) %>% unique()

prioritized_tbl_oi =
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  filter(ligand %in% filtered_ligands) %>% 
  select(niche, sender, receiver, ligand,  receptor, ligand_receptor, prioritization_score) %>% 
  distinct() %>% 
  inner_join(top_ligand_receptor_niche_df) %>% 
  group_by(ligand) %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(2, prioritization_score) %>% 
  ungroup() 
## colorの設定
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
  c("lavender") %>% 
  magrittr::set_names(
    prioritized_tbl_oi$receiver %>% 
      unique() %>% 
      sort())

circos_output = 
  make_circos_lr(
    prioritized_tbl_oi, 
    colors_sender, 
    colors_receiver)
