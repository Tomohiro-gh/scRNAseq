library(Seurat)
library(patchwork)
library(openxlsx)


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
########   Find all markers + Heatmap
### Description #####################
### listには使えない 'DefaultAssay<-' をクラス "list" のオブジェクトに適用できるようなメソッドがありません 
# Usage: fun.all.markers(seuratobject, assayname, samplename)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/SeuratFunctions.R
## version 22/5/18
fun.all.markers.v2 <- function(seuratobject, assayname,minimalpct, logfc ,samplename){
  
  DefaultAssay(seuratobject) <- assayname
  
  all.mk <- FindAllMarkers(seuratobject,
                           assay = assayname,
                           only.pos = TRUE, 
                           min.pct = minimalpct, 
                           logfc.threshold = logfc)
  # Findall markerで取り出したobjectを整理する
  TopAll = all.mk %>% 
    group_by(cluster) %>% 
    filter(p_val_adj < 0.01)
  
  # topallで，他のクラスタにないmarkerのみを抜き出す
  Unique = TopAll %>%
    filter(p_val_adj < 0.001) %>%
    group_by(gene) %>% 
    filter(n() == 1) %>% 
    group_by(cluster) %>%
    pull(gene)
  
  # 重複した遺伝子のみを抽出．
  Ambiguous = TopAll %>%
    filter(!(gene %in% Unique))  %>%
    pull(gene)
  
  # marker geneのUnique性をUniqueMarkerという新しいカラムを作って入れる
  TopAll$UniqueMarker[TopAll$gene %in% Unique] <- "Unique"
  TopAll$UniqueMarker[TopAll$gene %in% Ambiguous] <- "Ambiguous"
  
  # bookに保存
  filename <- paste0("FindAllMarkers_", samplename, ".xlsx")
  workbook <- createWorkbook()
  addWorksheet(workbook, sheetName = "All") # sheetを作る
  writeData(workbook, sheet = 1, x=TopAll, rowNames = FALSE) # dataを書き込む
  saveWorkbook(workbook, file = filename, overwrite = TRUE) # save
  
  ##### Heatmap section   
  # Gene top5を取り出してHeatmapを書いてみる si:とzgc:は除く
  require(stringr) #str_detectに必要
  
  Top5 <- TopAll %>% 
    group_by(cluster) %>% 
    filter(!str_detect(gene,"si:") & !str_detect(gene,"zgc:")) %>%
    slice_head(n = 5) %>% 
    pull(gene)
  
  pngname <- paste0("HeatmapTop5", samplename, ".png")
  DoHeatmap(object = seuratobject,
            features = Top5,
            slot = "scale.data") + 
    theme(axis.text.y = element_text(color = "black", size = 10))
  ggsave(pngname, dpi=300, width = 14, height = 10)
  
  # Unique gene top5を取り出してHeatmapを書いてみる si:とzgc:は除く
  Top5uni <- TopAll %>% 
    group_by(cluster) %>% 
    filter(UniqueMarker=="Unique") %>% 
    filter(!str_detect(gene,"si:") & !str_detect(gene,"zgc:")) %>%
    slice_head(n = 5) %>% 
    pull(gene)
  
  pngname <- paste0("HeatmapTop5Unique", samplename, ".png")
  DoHeatmap(object = seuratobject,
            features = Top5uni,
            slot = "scale.data") + 
    theme(axis.text.y = element_text(color = "black", size = 10))
  ggsave(pngname, dpi=300, width = 14, height = 10)
}

## Example (listをforを使ってやるのはできない)
# fun.all.markers(standard.un, "RNA", "standard.unwounded")
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


fun.all.markers <- function(seuratobject, assayname, samplename){
  
  DefaultAssay(seuratobject) <- assayname
  
  all.mk <- FindAllMarkers(seuratobject,
                           assay = assayname,
                           only.pos = TRUE, 
                           min.pct = 0.2, 
                           logfc.threshold = 0.25)
  # Findall markerで取り出したobjectを整理する
  TopAll = all.mk %>% 
    group_by(cluster) %>% 
    filter(p_val_adj < 0.01)
  
  # topallで，他のクラスタにないmarkerのみを抜き出す
  Unique = TopAll %>%
    filter(p_val_adj < 0.001) %>%
    group_by(gene) %>% 
    filter(n() == 1) %>% 
    group_by(cluster) %>%
    pull(gene)
  
  # 重複した遺伝子のみを抽出．
  Ambiguous = TopAll %>%
    filter(!(gene %in% Unique))  %>%
    pull(gene)
  
  # marker geneのUnique性をUniqueMarkerという新しいカラムを作って入れる
  TopAll$UniqueMarker[TopAll$gene %in% Unique] <- "Unique"
  TopAll$UniqueMarker[TopAll$gene %in% Ambiguous] <- "Ambiguous"
  
  # bookに保存
  filename <- paste0("FindAllMarkers_", samplename, ".xlsx")
  workbook <- createWorkbook()
  addWorksheet(workbook, sheetName = "All") # sheetを作る
  writeData(workbook, sheet = 1, x=TopAll, rowNames = FALSE) # dataを書き込む
  saveWorkbook(workbook, file = filename, overwrite = TRUE) # save
  
  ##### Heatmap section   
  # Gene top5を取り出してHeatmapを書いてみる si:とzgc:は除く
  require(stringr) #str_detectに必要
  
  Top5 <- TopAll %>% 
    group_by(cluster) %>% 
    filter(!str_detect(gene,"si:") & !str_detect(gene,"zgc:")) %>%
    slice_head(n = 5) %>% 
    pull(gene)
  
  pngname <- paste0("HeatmapTop5", samplename, ".png")
  DoHeatmap(object = seuratobject,
            features = Top5,
            slot = "scale.data") + 
    theme(axis.text.y = element_text(color = "black", size = 10))
  ggsave(pngname, dpi=300, width = 14, height = 10)
  
  # Unique gene top5を取り出してHeatmapを書いてみる si:とzgc:は除く
  Top5uni <- TopAll %>% 
    group_by(cluster) %>% 
    filter(UniqueMarker=="Unique") %>% 
    filter(!str_detect(gene,"si:") & !str_detect(gene,"zgc:")) %>%
    slice_head(n = 5) %>% 
    pull(gene)
  
  pngname <- paste0("HeatmapTop5Unique", samplename, ".png")
  DoHeatmap(object = seuratobject,
            features = Top5uni,
            slot = "scale.data") + 
    theme(axis.text.y = element_text(color = "black", size = 10))
  ggsave(pngname, dpi=300, width = 14, height = 10)
}
