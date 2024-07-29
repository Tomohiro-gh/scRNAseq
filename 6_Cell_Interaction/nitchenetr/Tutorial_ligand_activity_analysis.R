## ligandの発現を調べてみる

library(nichenetr)
library(tidyverse)

## Single-cell NicheNet’s ligand activity analysis

hnscc = readRDS("/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/hnscc_expression.rds")
class(hnscc)


expression = hnscc_expression$expression
sample_info = hnscc_expression$sample_info


tumors_remove = 
  c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids = 
  sample_info %>% 
  filter(`Lymph node` == 0) %>% 
  filter((tumor %in% tumors_remove == FALSE)) %>% 
  filter(`non-cancer cell type` == "CAF") %>% 
  .$cell

malignant_ids = 
  sample_info %>% 
  filter(`Lymph node` == 0) %>% 
  filter(`classified  as cancer cell` == 1) %>% 
  filter((tumor %in% tumors_remove == FALSE)) %>% 
  .$cell

expressed_genes_CAFs = 
  expression[CAF_ids,] %>% 
  apply(2,function(x){10*(2**x - 1)}) %>% 
  apply(2,function(x){log2(mean(x) + 1)}) %>% 
  .[. >= 4] %>% 
  names()

expressed_genes_malignant = 
  expression[malignant_ids,] %>% 
  apply(2,function(x){10*(2**x - 1)}) %>% 
  apply(2,function(x){log2(mean(x) + 1)}) %>% 
  .[. >= 4] %>% 
  names()



# Preparation:
location="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/"
#ligand_target_matrix = readRDS(paste0(location,"ligand_target_matrix.rds"))
lr_network = readRDS(paste0(location,"lr_network.rds"))
#weighted_networks= readRDS(paste0(location,"weighted_networks.rds"))


ligands = 
  lr_network$from %>% 
  unique()
expressed_ligands = 
  intersect(
    ligands,
    expressed_genes_CAFs)
receptors = 
  lr_network$to %>% 
  unique()
expressed_receptors = 
  intersect(
    receptors,
    expressed_genes_malignant)

potential_ligands = 
  lr_network %>% 
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  .$from %>% 
  unique()

head(potential_ligands)
#

# In a second step, we will scale the single-cell expression data (including only expressed genes).
background_expressed_genes = 
  expressed_genes_malignant %>% 
  .[. %in% rownames(ligand_target_matrix)]
expression_scaled = 
  expression %>% 
  .[malignant_ids,background_expressed_genes] %>% 
  scale_quantile()

## Ligand prioritization by regression analysis
normalized_ligand_activities = 
  normalize_single_cell_ligand_activities(ligand_activities)
# Furthermore, we will also show how you can perform additional analyses by linking the ligand activity in cells to other properties of cells in order to prioritize ligands. As toy example, we will score malignant cells here on the extent to which they express the core p-EMT gene “TGFBI”.
cell_scores_tbl = 
  tibble(
    cell = malignant_hn5_ids, 
    score = expression_scaled[malignant_hn5_ids,"TGFBI"])

#To do so, we frist need to process and normalize the ligand activities (i.e. pearson correlation values) to make different cells comparable. Here we use modified z-score normalization.
normalized_ligand_activities = 
  normalize_single_cell_ligand_activities(ligand_activities)


output_correlation_analysis = 
  single_ligand_activity_score_regression(
    normalized_ligand_activities,
    cell_scores_tbl)
output_correlation_analysis %>% 
  arrange(-pearson_regression) %>% 
  select(pearson_regression, ligand)

## # A tibble: 131 x 2
##    pearson_regression ligand  
##                 <dbl> <chr>   
##  1              0.525 TNC     
##  2              0.497 TFPI    
##  3              0.491 SEMA5A  
##  4              0.488 ANXA1   
##  5              0.473 TNFSF13B
##  6              0.462 IBSP    
##  7              0.449 HDGF    
##  8              0.443 HSP90B1 
##  9              0.431 CALM3   
## 10              0.428 CXCL11  
## # ... with 121 more rows


##. Visualize the relation between ligand activity and the cell's property score of interest
inner_join(
  cell_scores_tbl,
  normalized_ligand_activities) %>%
  ggplot(aes(score,TNC)) + 
  geom_point() + 
  geom_smooth(method = "lm")
