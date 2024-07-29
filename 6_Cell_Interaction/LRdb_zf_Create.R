# Creating Ligand-Receptor database in zebrafish

library(openxlsx)


# Ligand receptor interaction database
LR_human <- read.xlsx("/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/LRdb_SingleCellSignalR.xlsx", rowNames = FALSE, startRow = 1)


# Ortholog list
OrthologList <- read.xlsx("/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/ensembl_biomart/BiomaRt_ensembl_zf_hsa_ortholog_list_raw(41965).xlsx", rowNames = FALSE, startRow = 1)　
colnames(OrthologList) <- c("zf_emsebl","zf_symbol", "human_ensembl","human_symbol")
  dim(OrthologList) # 41965     4


for(i in 1:nrow(LR_human)){
  # i 行目を取り出し
  r <- LR_human[,i]
  
  # 変換　human to zebrafish
  gene.zf <- OrthologList %>%
    filter(human_symbol %in% r)  %>%
      pull(zf_symbol)
  
  # 2つ以上であれば
  # dataframeへ追加
  
  
}