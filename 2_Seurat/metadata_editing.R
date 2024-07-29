## Seurat  meta dataの編集

## ２つのidentを繋ぎ合わせた文字列で新たにidentを作成する
SeuratObject$NEWmetadata <- paste(Idents(SeuratObject), SeuratObject$metadata1, sep = "_")
Idents(SeuratObject) <- "NEWmetadata"
as.data.frame(table(Idents(SeuratObject)))
## 数を確認できる

## 文字を削ってmetadataへ入れる
substring(SeuratObject$metadata1, 4,7)  #1から2文字目までを取り出す
Idents(SeuratObject) <- "ClusterxStage"
as.data.frame(table(Idents(SeuratObject)))



## >> 06/13/23 
## metadataの一部を編集する．
## 例えば，MC をPericyteとSMCに分ける
## Original codes: 230612_ClusterAnnottion_Dim26Res0.5.R
#"/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp155_reanalysis_of_2022DC_LCangiogenesis/3_ClusterAnnotation_Dim26Res0.5/230612_ClusterAnnottion_Dim26Res0.5.R" 


## Subclusterのところだけannotateionを変える
seurat_obj@meta.data <- seurat_obj@meta.data  %>% 
  mutate(CellAnnotation_Fine = CellAnnotation) %>%
  mutate(CellAnnotation_Fine = case_when(
    rownames(.) %in% Barcodes_EC ~ "Endothelial Cell",
    rownames(.) %in% Barcodes_LEC ~ "Lymphatic EC",
    rownames(.) %in% Barcodes_PC ~ "Pericyte",
    rownames(.) %in% Barcodes_SMC ~ "SMC",
    rownames(.) %in% Barcodes_Epi_unknown ~ "EpiECMC?",
    TRUE ~ .$CellAnnotation))