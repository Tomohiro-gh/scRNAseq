library(SingleCellSignalR)
library(Seurat)


wd = ("/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp170_LungEC_Aging_scRNAseq/CellCellInteraction/2_liana/230516_SingleCellSignalR_from_liana_results")
setwd(wd)

## 参考になるデータ
cell_signal = readRDS("/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp170_LungEC_Aging_scRNAseq/CellCellInteraction/3_SingleCellSignalR/cell_signaling_230513.rds")
colnames(cell_signal$`12w-12w`)
# "12w"              "12w"              "interaction type" "LRscore"         
## Lianaの結果のLR スコアから，SingleCellSignalRのCircus plotが書けないか？


## loading liana resutls
liana_res = readRDS("/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp170_LungEC_Aging_scRNAseq/CellCellInteraction/2_liana/liana_res_230513.rds")

## data frameの整形
# col: ligand, receptor, interaction type, LRscoreにまとめる

liana_res
colnames(liana_res)
# [1] "source"                 "target"                 "c"        
# [4] "receptor.complex"       "aggregate_rank"         "mean_rank"             
# [7] "natmi.edge_specificity" "natmi.rank"             "connectome.weight_sc"  
# [10] "connectome.rank"        "logfc.logfc_comb"       "logfc.rank"            
# [13] "sca.LRscore"            "sca.rank"               "cellphonedb.pvalue"    
# [16] "cellphonedb.rank"       "Ligand_Receptor"       

## ligand complex, receptor complex and 



liana_res_yy <- 
  liana_res %>%
  filter(source==c("12w") & target == c("12w")) 

liana_res_oo <- 
  liana_res %>%
  filter(source==c("72w") & target == c("72w")) 


## single cell signalRは SCAに格納
## 13の LR scoreを使えばいい
liana_res_yy_top30 <- 
  liana_res %>%
  filter(source==c("12w") & target == c("12w")) %>%
  slice_min(aggregate_rank, n=30) %>% 
  select(c("ligand.complex", "receptor.complex", "sca.LRscore")) %>%
  mutate("interaction type" = "autocrine") %>% 
  select(ligand.complex, receptor.complex, "interaction type", sca.LRscore)
colnames(liana_res_yy_top30) <- c("12w", "12w", "interaction\ type", "LRscore")

## 13の LR scoreを使えばいい
liana_res_oo_top30 <- 
  liana_res %>%
  filter(source==c("72w") & target == c("72w")) %>%
  slice_min(aggregate_rank, n=30) %>% 
  select(c("ligand.complex", "receptor.complex", "sca.LRscore")) %>%
  mutate("interaction type" = "autocrine") %>% 
  select(ligand.complex, receptor.complex, "interaction type", sca.LRscore)
colnames(liana_res_oo_top30) <- c("72w", "72w", "interaction\ type", "LRscore")


cell_signal_from_liana = list(
  `12w-12w` = as.data.frame(liana_res_yy_top30),
  `72w-72w` = as.data.frame(liana_res_oo_top30)
)

cell_signal_from_liana

names(cell_signal_from_liana)

## SinglecellSignalRを動かしてみる
## 
visualize_interactions(signal = cell_signal_from_liana, write.in = c(1), limit=30)
visualize_interactions(signal = cell_signal_from_liana, write.in = c(2), limit=30)
visualize_interactions(signal = cell_signal_from_liana, write.in = c(1), limit=20)
visualize_interactions(signal = cell_signal_from_liana, write.in = c(2), limit=20)

## 完成




## 次は oldもしくはyoungのみでみられたinteractionを書いてみる
# liana_res = readRDS("liana_res_230513.rds")

## ligand_receptor項目追加
liana_res <- 
  liana_res %>% 
  mutate(
    Ligand_Receptor = paste(ligand.complex,receptor.complex, sep=" -> ")
  )

liana_res_yy <- 
  liana_res %>%
  filter(source==c("12w") & target == c("12w")) 
dim(liana_res_yy) #341  17

liana_res_oo <- 
  liana_res %>%
  filter(source==c("72w") & target == c("72w")) 
dim(liana_res_oo) #312  17


common_int <- 
  intersect(liana_res_yy$Ligand_Receptor, liana_res_oo$Ligand_Receptor)
length(common_int) #296


## commonなやつは除く
liana_res_only_yy_top30 <- 
  liana_res %>%
  filter(source==c("12w") & target == c("12w")) %>%
  filter(!Ligand_Receptor %in% common_int) %>% 
  slice_min(aggregate_rank, n=30) %>% 
  select(c("ligand.complex", "receptor.complex", "sca.LRscore")) %>%
  mutate("interaction type" = "autocrine") %>% 
  select(ligand.complex, receptor.complex, "interaction type", sca.LRscore)
colnames(liana_res_only_yy_top30) <- c("12w", "12w", "interaction\ type", "LRscore")

## 13の LR scoreを使えばいい
liana_res_only_oo <- 
  liana_res %>%
  filter(source==c("72w") & target == c("72w")) %>%
  filter(!Ligand_Receptor %in% common_int) %>% 
  slice_min(aggregate_rank, n=30) %>% 
  select(c("ligand.complex", "receptor.complex", "sca.LRscore")) %>%
  mutate("interaction type" = "autocrine") %>% 
  select(ligand.complex, receptor.complex, "interaction type", sca.LRscore)
colnames(liana_res_only_oo) <- c("72w", "72w", "interaction\ type", "LRscore")

## For SingleCellSignalR
cell_signal_from_liana_only_yy_or_oo = list(
  `12w-12w` = as.data.frame(liana_res_only_yy_top30),
  `72w-72w` = as.data.frame(liana_res_only_oo)
)

visualize_interactions(signal = cell_signal_from_liana_only_yy_or_oo, write.in = c(1), limit=30)
visualize_interactions(signal = cell_signal_from_liana_only_yy_or_oo, write.in = c(2), limit=30)
visualize_interactions(signal = cell_signal_from_liana_only_yy_or_oo, write.in = c(1), limit=20)
# visualize_interactions(signal = cell_signal_from_liana_only_yy_or_oo, write.in = c(2), limit=20)
