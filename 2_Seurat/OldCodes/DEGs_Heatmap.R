library(Seurat)
library(sctransform)
library(glmGamPoi)
library(ggplot2)
library(SeuratData)
library(patchwork)
library(dplyr)
library(Matrix)
library(gdata)
library(reshape)
library(reshape2)
library(openxlsx)
library(tidyverse)


############# 05/23- DEG analysis across conditions ###############
wd = "/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp130_WH_scRNAseq_Analysis/SeuratAnalysis_v5/2-4_Integration_SCT_thenRBCremoval"
setwd(wd)


### Cluster annnotationとRBCを除いて解析
SCT22375 <- readRDS("22375Cells_SCT_Annotated.rds")
# Annotation済み，PrepSCTfindMarkers済み
# Confirmation
table(Idents(SCT22375))


#confirm cluster idents (modes)
levels(SCT22375) # 現在のactiveになっているmetadata
# "RBC"  "Myeloid" "Mesenchyme" "Epithelia"  "Osteoblast" "Neutrophil" "DC,Tcell"  
colnames(SCT22375@meta.data)

# ClusterとClusterxStageという２つの新しい identを作成する
SCT22375$Clusters <- Idents(SCT22375)
SCT22375$ClusterxStage <- paste(Idents(SCT22375), 
                                SCT22375$Condition, 
                                sep = "_")

### Conrolと各ステージを比べて DEGを算出する  
# IdentをClusterxStageにセットする
Idents(SCT22375) <- "ClusterxStage"

# FindMarkersで 現在の2つのクラスを比べ，変動が見られる遺伝子を見つける
# level.list <- levels(SCT22375$ClusterxStage) # ClusterxStage
fun.DEgenes <- function(object, cluster, name){
  
  stage <- c("Unwounded", "2dpi", "4dpi", "7dpi")
  
  Idents(object) <- "ClusterxStage"
  # Cluster x stageでのobjectをそれぞれ作る
  group <- NULL #最初に作っておかないとエラーが出る
  for(i in 1:length(stage)){
    group[i] <- paste(cluster, stage[i], sep="_")
  }
  
  workbook <- createWorkbook()
  
  # Ident.2にcotrolをおく logFCはident.2に対する発現の割合
  for(j in 1:(length(group)-1)){
    DE <- FindMarkers(object,
                      ident.1 = group[j+1],
                      ident.2 = group[1],
                      verbose = FALSE)
    DEgenes = DE %>% filter(p_val_adj < 0.05)
    
    #シートへ保存：追加OK    
    sName = paste(stage[1],"vs", stage[j+1], sep="")
    addWorksheet(workbook, sheetName = sName)
    writeData(workbook, sheet = sName, x=DEgenes,rowNames = TRUE)
    saveWorkbook(workbook, 
                 file = paste0("DEgenes_", cluster, "_", name, ".xlsx"),
                 overwrite = T)
  }
}

# Execution
#今回解析するクラスターたちを指定
cluster.list <- c("Myeloid", "Mesenchyme", "Epithelia","Osteoblast","Neutrophil", "DC,Tcell")
# すべてのクラスターで一度に行う．
for(a_cluster in cluster.list){
  fun.DEgenes(SCT22375, a_cluster, "22375Cells_SCT")
}

## OK !!! 
# Excel file が生成される "DEgenes_xxxxxx_22375Cells_SCT.xlsx" 
# これをinputにしてheatmapを描く



## SCale dataに全ての遺伝子を入れておく
DefaultAssay(SCT22375) <- "SCT"
SCT22375 <- ScaleData(SCT22375, features = rownames(SCT22375))
nrow(SCT22375[["SCT"]]@scale.data) #22018 OK


## Angiogenesis
## 上で出したexel fileを元に，DEGをHeatmapで表示する
# AngiogenesisのGOにannotateされている遺伝子を読み込む
GO0001525 <- read.xlsx("/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/GO_Genes_zf/GO0001525_zf_Angiogenesis.xlsx",rowNames = FALSE, startRow = 2)　
head(GO0001525)
genes_Angiogenesis <- GO0001525$Symbol
length(genes_Angiogenesis) #373
# genes_Angiogenesisへgeneを格納（symbol）


#####  Excel fileへ読込む excel fileを格納する．
files <- NULL
for (i in 1:length(cluster.list)){
  files[i] <- paste0("DEgenes_", cluster.list[i], "_22375Cells_SCT.xlsx")
}


###### クラスターの中で発現変動遺伝子をheatmapで表示する
fun.DEGs.HP <- function(obj, datapath, cluster){
  sheets=c(1,2,3)
  stage=c("vs2dpi","vs4dpi","vs7dpi")
  
  workbook <- createWorkbook()
  
  for (k in 1:length(sheets)){
    exceldata <- read.xlsx(datapath, sheet=sheets[k], rowNames = FALSE)
    # colnamesの変更はdataframe全体
    ncol(exceldata)
    colnames(exceldata) <- c("Symbol", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
    DEGlist <- exceldata %>%
      filter(p_val_adj < 0.05) %>% 
      filter(Symbol %in% genes_Angiogenesis)
    
    addWorksheet(workbook, sheetName = stage[k])
    writeData(workbook, sheet = sheets[k], x=DEGlist, rowNames = FALSE)
    if(k==1){
      df <- as.data.frame(DEGlist)}
    else{df <- rbind(df, DEGlist)}
    DEGlist <- df %>%
      filter(! duplicated(Symbol))
  } 
  
  file_name <- paste0("DEGenesAngiogenesis_in_", cluster,".xlsx")
  saveWorkbook(workbook, file = file_name, overwrite = TRUE)
  
  DEGs <- DEGlist[ ,1]
  #subset関数でクラスターを抽出: Heatmapを書くためには，RNAassayでscale dataしている必要あり
  if(length(DEGs)>=1){
    a_cluster <- subset(x = obj, idents = cluster)
    Idents(a_cluster) <- "orig.ident"
    
    #Heatmapの作成
    file_name2 <- paste0("DEGenesAngiogenesis_in_", cluster,"_Heatmap.png")
    png(file_name2, w=15, h=10, units="in", res=200)
    dh <- DoHeatmap(object = a_cluster,
                    features = DEGs) + 
      theme(axis.text.y = element_text(color = "black", size = 10))
    plot(dh)
    dev.off()
  }
  
}

# execute one file (Myeloid)
fun.DEGs.HP(wh.int.woRBC.61.15, files[13], Clusters[13])

# execution all 
Idents(SCT22375) <- "Clusters"
for (i in 1:length(files)){
  fun.DEGs.HP(SCT22375, files[i], cluster.list[i])
}
#### Successfully run !


