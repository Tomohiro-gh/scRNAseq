# Seurat function and packages

library(Seurat)
library(SeuratData) # for integration
library(ggplot2)
library(patchwork)
library(cowplot)
# use below packages in one sample analysis
library(dplyr)
library(Matrix)
library(gdata)
library(reshape2)
library(openxlsx)
library(tidyverse)
library(stringr)

####### Standard Normalization :
## https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
fun.normalize.PCA <- function(seuratobject){
  # scale factor 10,000でnormalization
  seuratobject <- NormalizeData(seuratobject, 
                                normalization.method = "LogNormalize", 
                                scale.factor = 10000)
  # 変動遺伝子を見つける
  seuratobject <- FindVariableFeatures(seuratobject, 
                                       selection.method = "vst", 
                                       nfeatures = 2000)
  # a linear transformation (‘scaling’) that is a standard pre-processing step 
  # prior to dimensional reduction techniques like PCA.
  all.genes <- rownames(seuratobject)
  seuratobject <- ScaleData(seuratobject, features=all.genes)
  
  seuratobject <- RunPCA(seuratobject, 
                         features = VariableFeatures(object = seuratobject),
                         verbose = FALSE)
  
  return(seuratobject)
  
}


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
######### Dim の決定 : elbow plot
fun.elb.plot <- function(object,filename){
  elb <- ElbowPlot(object, ndims = 50) + 
    labs(title=paste0("ElbowPlot : ", filename)) + 
    theme(plot.title = element_text(hjust = 0.5))
  plot(elb)
  
  ggsave(paste0("PC_ElbowPlot_", filename,".png"), width=7, height=5, dpi=300)
}
# Execution
fun.elb.plot(woRBC8898_int, "1355Cells")
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




# >>>>>>>>>>> Cluster finding >>>>>>>>>>>>>>>>>>>>> 05/23
# This function is requiered for predefined dim, res range like above
# PCAまで行っているものが対象
# Usage: fun.cluster.finding.UMAP(sub, filename)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/SeuratFunctions.R
con_dim = c(26,31,38,46) # variable depending on condition
con_res = seq(0.1, 0.3, by=0.05) # variable depending on condition
fun.cluster.finding.UMAP <- function(sub, filename){
  pdfname <- paste0("ClusteringVariationUMAP_", filename, ".pdf")
  
  pdf(pdfname, width = 6, height = 6)
  
  for (i in 1:length(con_dim)){
    sub <- RunUMAP(sub, dims = 1:con_dim[i]) #3
    sub <- FindNeighbors(object = sub, dims = 1:con_dim[i]) #4
    for (j in 1:length(con_res)){
      sub <- FindClusters(object = sub, resolution = con_res[j]) #5
      u <- UMAPPlot(object =sub , label=T) + 
        NoLegend() + 
        labs(title=paste0("Dimention = ", con_dim[i], ",  Resolution = ", con_res[j])) + 
        theme(plot.title = element_text(hjust = 0.5))
      plot(u) 
    }
  }
  dev.off()
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





### 
# Clusteringの条件が決まったら，オブジェクトを作るために実行
fun.pca.tsne.umap <- function(object, dim, res){
  DefaultAssay(object) <- "integrated"
  object2 <- object
  object2 <- RunPCA(object2, verbose = FALSE)
  object2 <- RunUMAP(object2, reduction = "pca", dims = 1:dim)
  object2 <- RunTSNE(object2, reduction = "pca", dims = 1:dim)
  object2 <- FindNeighbors(object2, reduction = "pca", dims = 1:dim)
  object2 <- FindClusters(object2, resolution = res)
  
  return(object2)
}
#Execution (新たなobjectをつくる（noramlizedした後）
#wh.int.38.15 <- fun.clustering(wh.int, 38, 0.15)



######### Clustering 条件検討 2/2
# Clusteringに使用するDimとresolutionの検討
fun.det.dim <- function(sub, filename,dim){
    pdfname <- paste("ClusteringCondition(1of2)-", filename, ".pdf", sep = "")
    sub2 <- sub
  
  DefaultAssay(sub2) <- "integrated"
  sub2 <- ScaleData(sub2, verbose=FALSE) #1
  sub2 <- RunPCA(sub2, npcs=dim, verbose=FALSE) #2
  
  # 1) elbow plot : 降下曲線から水平になる点（Elbow）が最適な主成分の次元数
  e <- ElbowPlot(sub2, ndims = dim) + 
    labs(title=paste("elbow plot : ", filename, sep="")) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  # 2)  heatmap : 二分されているように見えなければ、その主成分にはもう情報がないと言える
  pdf(pdfname, width = 10, height = 8)
  plot(e)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
  DimHeatmap(sub2, dims = 1:15, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = filename, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
  DimHeatmap(sub2, dims = 16:30, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = filename, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
  DimHeatmap(sub2, dims = 31:45, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = filename, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
  DimHeatmap(sub2, dims = 46:60, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = filename, cex=1)
  
  dev.off()
}
# Execution
# fun.det.dim(wh.int, "Integration",60)


######### Clustering 条件検討 2/2
# Clusteringに使用するDimとresolutionの検討
fun.cluster.condition <- function(sub, filename){
  pdfname <- paste("ClusteringCondition(2of2)-", filename, ".pdf", sep="")
  sub2 <- sub
  con_dim = c(26,31,38,46) # variable depending on condition
  con_res = seq(0.1, 0.3, by=0.05) # variable depending on condition
  
  DefaultAssay(sub2) <- "integrated"
  sub2 <- ScaleData(sub2, verbose=FALSE) #1
  sub2 <- RunPCA(sub2, npcs=50, verbose=FALSE) #2
  
  pdf(pdfname, width = 6, height = 6)
  
  for (i in 1:length(con_dim)){
    sub2 <- RunUMAP(sub2, dims = 1:con_dim[i]) #3
    sub2 <- FindNeighbors(object = sub2, dims = 1:con_dim[i]) #4
    for (j in 1:length(con_res)){
      sub2 <- FindClusters(object = sub2, resolution = con_res[j]) #5
      u <- UMAPPlot(object =sub2 , label=T) + 
        NoLegend() + 
        labs(title=paste("Dimention = ", con_dim[i], ",  Resolution = ", con_res[j]), sep="") + 
        theme(plot.title = element_text(hjust = 0.5))
      plot(u) 
    }
  }
  dev.off()
}
#Execcution
#fun.cluster.condition(wh.int, "Integration")





# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
########   Find all markers + Heatmap
### Description #####################
### listには使えない 'DefaultAssay<-' をクラス "list" のオブジェクトに適用できるようなメソッドがありません 
# Usage: fun.all.markers(seuratobject, assayname, samplename)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/SeuratFunctions.R
## version 22/5/18
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
    filter(p_val_adj < 0.000001) %>%
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
fun.all.markers(standard.un, "RNA", "standard.unwounded")
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## DEgenesを同定する = 条件で変動する遺伝子をFindMarkerで見つける
# Usage: fun.DEgenes(seuratobject, assayname, samplename)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/SeuratFunctions.R
## version 22/5/18

## preparation : # Cluster(番号 or name)  x Stageの identを作る
  Idents(object) <- "seurat_clusters"
  msen20.020$ClusterxStage <- paste(Idents(object), 
                                  object$Condition, 
                                  sep = "_")

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


# Example : １度に全てのclusterについて解析を行いたい場合．
Idents(object) <- "seurat_clusters"
cluster.list　<- levels(object)

for(a_cluster in cluster.list){
  fun.DEgenes(object, a_cluster, "Mesenchyme1789Cells_SCT_Dim20Res0.2")
}


# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
