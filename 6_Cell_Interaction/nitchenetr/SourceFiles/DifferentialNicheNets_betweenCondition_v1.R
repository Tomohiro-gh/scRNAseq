## Differential NicheNet analysis between conditions of interest
## https://github.com/saeyslab/nichenetr/blob/master/vignettes/differential_nichenet_pEMT.md

# ２つのcondition間でのDifferential NicheNetやるためのまとめcode
# 実行方法は，
# source('/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/SourceFiles/DifferentialNicheNets_betweenCondition.R')

library(nichenetr)
library(RColorBrewer)
library(tidyverse)
library(Seurat) #
# additional libraries
library(stringr)
library(dplyr)
library(RColorBrewer)

##定義しておくべきもの #1 
## Seuratobject -> 必ずseurat_objで入れておく
seurat_obj 


## 定義しておくべきもの #2 
## Idents -> cluster (celltype x conditionなど)のindentを作り，セットしておく
## example )
# seurat_obj@meta.data$State_Annotation1 = 
#.  paste(seurat_obj@meta.data$State, 
#      seurat_obj@meta.data$Annotation1,
#      sep = "_") 

# celltype_id = "State_Annotation1" # metadata column name of the cell type of interest
#.  seurat_obj = 
#  SetIdent(
#    seurat_obj,
#    value = seurat_obj[[celltype_id]])


## 定義しておくべきもの #3
## niche list: 一番大事

#niches = list(
 # To_EC_in_UW = list(
  #  "sender" = c("UW_PC",
   #              "UW_PC",
    #             "UW_MC int.",
  #               "UW_SMC"),
  #  "receiver" = c("UW_EC")),
  #To_EC_in_WO = list(
  #  "sender" = c("WO_PC",
  #               "WO_SMC",
  #               "WO_MC int.",
  #               "WO_Mitotic EC",
  #               "WO_Mitotic MC"),
  #  "receiver" = c("WO_EC"))
#) # user adaptation required on own dataset


## 定義しておくべきもの #4
# organism = "mouse"
# or
# organism = "human"


## Pre defined parameter:
  #1) expression_pct = 0.10 (STEP2)
  #2) lfc_cutoff = 0.15 (STEP4)　# recommended for 10x as min_lfc cutoff. 
  #3) top_n_target = 250 (STEP4)
 

### Read in the NicheNet ligand-receptor network and ligand-target matrix
location="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/"
ligand_target_matrix = 
  readRDS(paste0(location,"ligand_target_matrix.rds"))

# ligand receptor network
lr_network = readRDS(paste0(location,"lr_network.rds"))
lr_network = 
  lr_network %>% 
  mutate(bonafide = ! database %in% 
           c("ppi_prediction","ppi_prediction_go"))
lr_network = 
  lr_network %>% 
  dplyr::rename(ligand = from, receptor = to) %>% 
  distinct(ligand, receptor, bonafide)

# organism = "mouse" # or human
# mouse遺伝子へ変換
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

# dim(lr_network)　# A tibble: 11592     3
# dim(ligand_target_matrix) # 17330   644

## 1 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Define the niches/microenvironments of interest

##  2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate differential expression between the niches
assay_oi = "RNA" # other possibilities: RNA,...

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

print("STEP2 finished")
## ３ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##  Optional: Calculate differential expression between the different spatial regions
include_spatial_info_sender = FALSE 
include_spatial_info_receiver = FALSE #
# if not spatial info to include: put this to false # user adaptation required on own dataset

include_spatial_info_sender = FALSE # if not spatial info to include: put this to false # user adaptation required on own dataset
include_spatial_info_receiver = FALSE # if spatial info to include: put this to true # user adaptation required on own dataset
## ここは　skip

#spatial_info = 
# tibble(
#    celltype_region_oi = "CAF_High", 
#    celltype_other_region = "myofibroblast_High", 
#    niche =  "pEMT_High_niche", 
#    celltype_type = "sender") # user adaptation required on own dataset

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

print("STEP3 finished")


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



## それぞれのnicheでのtarget geneを算出
geneset_niche <- list()

for (i in 1:length(niches)){
  a_geneset_niche = 
    DE_receiver_processed_targets %>% 
    filter(
      receiver == niches[[i]]$receiver & 
        target_score >= lfc_cutoff & 
        target_significant == 1 & 
        target_present == 1) %>% 
    pull(target) %>% 
    unique()
  
  # listの結合
  geneset_niche <- c(geneset_niche, list(a_geneset_niche))
  
  geneset_niche[i] %>% 
    setdiff(rownames(ligand_target_matrix)) %>% 
    print()
  
}

print(dplyr::glimpse(geneset_niche))

print("finish geneset_niche list")

## 数的には少ないほうなのでこれで行く
top_n_target = 250

niche_geneset_list <- list()
niche_geneset_list_name <- NULL

for(i in 1:length(niches)){
  
  a_geneset_niche = 
    geneset_niche[i] %>% 
    unlist()
  
  a_niche_geneset_list = list(
    "receiver" = niches[[i]]$receiver,
    "geneset" = a_geneset_niche,
    "background" = background)

  niche_geneset_list <-
    append(niche_geneset_list,
           list(a_niche_geneset_list))
  
  niche_geneset_list_name[i] <- names(niches[i])
  
  }

names(niche_geneset_list) <- niche_geneset_list_name
 print("check niche_geneset_list")
 dplyr::glimpse(niche_geneset_list)
 
 
## 元のtrailでは以下のように作る 
#  list1 = list(
#    "receiver" = niches[[1]]$receiver,
#    "geneset" = geneset_niche[1],
#    "background" = background),
#  list2 = list(
#    "receiver" = niches[[2]]$receiver,
#    "geneset" = geneset_niche[2],
#    "background" = background),
#  list3 = list(
#    "receiver" = niches[[3]]$receiver,
#    "geneset" = geneset_niche[3],
#    "background" = background),
#  list4 = list(
#    "receiver" = niches[[4]]$receiver,
#    "geneset" = geneset_niche[4],
#    "background" = background)
#)

ligand_activities_targets = 
  get_ligand_activities_targets(
    niche_geneset_list = niche_geneset_list, 
    ligand_target_matrix = ligand_target_matrix, 
    top_n_target = top_n_target)

print("STEP4 finished")


## 5 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Calculate (scaled) expression of ligands, receptors and targets
features_oi = 
  union(lr_network$ligand, lr_network$receptor) %>% 
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
  rename(celltype = id, 
         gene = features.plot, 
         expression = avg.exp, 
         expression_scaled = avg.exp.scaled, 
         fraction = pct.exp) %>%
  mutate(fraction = fraction/100) %>% 
  as_tibble() %>% 
  select(celltype, 
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
    scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% 
  mutate(
    ligand_fraction_adapted = ligand_fraction) %>% 
  mutate_cond(
    ligand_fraction >= expression_pct,
    ligand_fraction_adapted = expression_pct)  %>% 
  mutate(
    scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor = 
  exprs_tbl_receptor %>% 
  mutate(
    scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% 
  mutate(
    receptor_fraction_adapted = receptor_fraction) %>% 
  mutate_cond(
    receptor_fraction >= expression_pct, 
    receptor_fraction_adapted = expression_pct)  %>% 
  mutate(
    scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))


print("STEP5 finished")

## 6 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Expression fraction and receptor 
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

print("STEP6 finished")

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

## Outputはリストにして
output = list(
  DE_sender_receiver = DE_sender_receiver, 
  ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
  sender_spatial_DE_processed = sender_spatial_DE_processed, 
  receiver_spatial_DE_processed = receiver_spatial_DE_processed,
  ligand_activities_targets = ligand_activities_targets, 
  DE_receiver_processed_targets = DE_receiver_processed_targets, 
  exprs_tbl_ligand = exprs_tbl_ligand,  
  exprs_tbl_receptor = exprs_tbl_receptor, 
  exprs_tbl_target = exprs_tbl_target)

# prioritization_tables
prioritization_tables = 
  get_prioritization_tables(output, prioritizing_weights)


## Confirmation
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


#outputとprioritization_tablesをsave
saveRDS(output, "1_output.rds")
saveRDS(prioritization_tables, "2_prioritization_tables.rds")

print("STEP7 finished")


## 8 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Visualization of the Differential NicheNet output
top_ligand_niche_df = 
  prioritization_tables$prioritization_tbl_ligand_receptor %>% 
  select(niche, sender, receiver, ligand, receptor, prioritization_score) %>% 
  group_by(ligand) %>% 
  top_n(1, prioritization_score) %>% 
  ungroup() %>% 
  select(ligand,
         receptor,
         niche) %>% 
  rename(top_niche = niche)

top_ligand_receptor_niche_df =
  prioritization_tables$prioritization_tbl_ligand_receptor %>%
  select(niche, sender, receiver, ligand, receptor,prioritization_score) %>%
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



## FUNCTION to output
graphics.output <- function(filename, my_receiver_oi){

receiver_oi = my_receiver_oi

filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  pull(ligand) %>% 
  unique()

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

lfc_plot = 
  make_ligand_receptor_lfc_plot(
    receiver_oi, prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor,
    plot_legend = FALSE,
    heights = NULL,
    widths = NULL)

pdf(paste0(filename, "_lfc_plot.pdf"),
    width=14, height=20)
plot(lfc_plot)
dev.off()


# Ligand expression, activity and target genes
# Active target gene inference - cf Default NicheNet
#Now: visualization of ligand activity and ligand-target links
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

## 
pdf(paste0(filename, "_exprs_activity_target_plot.pdf"),
    width=20, height=10)
plot(exprs_activity_target_plot$combined_plot)
plot(exprs_activity_target_plot$legends)

dev.off()


# If this plot contains too much information because we look at many hits (top 50 ligands), you can make this plot of course for less ligands as well, eg for the top20.
filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>%
  top_n(20, prioritization_score) %>% 
  pull(ligand) %>% 
  unique()

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

exprs_activity_target_plot = 
  make_ligand_activity_target_exprs_plot(
    receiver_oi, prioritized_tbl_oi, 
    prioritization_tables$prioritization_tbl_ligand_receptor, 
    prioritization_tables$prioritization_tbl_ligand_target,
    output$exprs_tbl_ligand, 
    output$exprs_tbl_target,
    lfc_cutoff,
    ligand_target_matrix, 
    plot_legend = FALSE,
    heights = NULL,
    widths = NULL)

## 
pdf(paste0(filename, "_exprs_activity_target_plot_selected.pdf"),
    width=15, height=7)
plot(exprs_activity_target_plot$combined_plot)
dev.off()

## Circos plot of prioritized ligand-receptor pairs
# Because a top50 is too much to visualize in a circos plot, we will only visualize the top 15.
filtered_ligands = 
  ligand_prioritized_tbl_oi %>% 
  filter(receiver == receiver_oi) %>% 
  top_n(20, prioritization_score) %>% 
  pull(ligand) %>% 
  unique()

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

colors_sender =
  brewer.pal(
    n = prioritized_tbl_oi$sender %>%
      unique() %>%
      sort() %>% 
      length(), name = 'Accent') %>%
  magrittr::set_names(
    prioritized_tbl_oi$sender %>%
      unique() %>%
      sort())
colors_receiver = 
  c("tomato")  %>%
  magrittr::set_names(
    prioritized_tbl_oi$receiver %>%
      unique() %>%
      sort())

# save
pdf(paste0(filename, "_Circosplot.pdf"),
    width=6, height=6)

circos_output = 
  make_circos_lr(
    prioritized_tbl_oi, 
    colors_sender,
    colors_receiver)

dev.off()


pdf(paste0(filename, "_Circosplot_top.pdf"),
    width=6, height=6)

circos_output = 
  make_circos_lr(
    prioritized_tbl_oi, 
    colors_sender,
    colors_receiver)

dev.off()


}

print("STEP8 (Graphical output) finished")

##########
## 全てのpatternで保存する
wd_current = getwd()

for (i in 1:length(niches)){
  # nicheの名前
  a_name <- names(niches[i])
  # nicheの名前のディレクトリ．ここへ保存する
  dir.create(a_name)
  setwd(paste0(wd_current, "/", a_name))
  
  receiver_oi = niches[[i]]$receiver
  
  graphics.output(a_name, receiver_oi)
  
  rm(a_name)
  setwd(wd_current)
  
}


print("Done")