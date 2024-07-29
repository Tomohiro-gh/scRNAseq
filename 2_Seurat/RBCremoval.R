library(Seurat)
library(SeuratData)

#### RBCや特定の細胞のbarcodeだけをのこす，除く


#### RBCを取りたい場合　########
# Cluter を確認
table(Idents(wh.sct))

#ここからCluster 0,8 ("RBC"と名付けたもの)を取り出す
# or RBCと名前がつけられたものだけを取り出す
RBCell <- subset(x=wh.sct, idents=c("0","1"))

#RBCとつけられた名前だけを取り出す
RBC_name <- rownames(RBCell[[]])
saveRDS(RBC_name, "RBC_CellBarcode_SCT.rds")
#RBC_name <- readRDS("RBC_CellBarcode.rds")

# whはmergeしただけのサンプルがいい
#RBC_nameに含まれる細胞を全て取り除き，新たなobjectを作成する
wh.woRBC <- wh[,!colnames(wh) %in% RBC_name]

wh.woRBC #confirmation
# An object of class Seurat 
# 20884 features across 6947 samples within 1 assay 
# Active assay: RNA (20884 features, 0 variable features)

saveRDS(wh.woRBC, "WHwoRBC_notNorm2.rds")



## 5/21 追記
#identsに -c をつけることで，指定したものだけを除くこともできる
stan.un.woRBC <- subset(x=stan.un, idents=-c("0","1"))

