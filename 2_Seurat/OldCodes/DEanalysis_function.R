

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
### One cluster -> DEG analysis across condition 
## 06/02/22
# Usage: fun.DEgenes.subsetdata(object, name)
## /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/DEanalysis_function.R
fun.DEgenes.subsetdata <- function(object, name){
  DefaultAssay(object) <- "RNA" # DEの時はRNAで
  
  stage <- c("Unwounded", "2dpi", "4dpi", "7dpi")
  # Cluster x stageでのobjectをそれぞれ作る
  #group <- NULL #最初に作っておかないとエラーが出る
  #for(i in 1:length(stage)){
  # group[i] <- paste(cluster, stage[i], sep="_")
  #}
  workbook <- createWorkbook()
  
  # Ident.2にcotrolをおく logFCはident.2に対する発現の割合
  for(j in 1:(length(stage)-1)){
    DE <- FindMarkers(object,
                      ident.1 = stage[j+1],
                      ident.2 = stage[1],
                      verbose = FALSE)
    DEgenes = DE %>% filter(p_val_adj < 0.05)
    
    #シートへ保存：追加OK    
    sName = paste0(stage[1], "vs", stage[j+1])
    addWorksheet(workbook, sheetName = sName)
    writeData(workbook, sheet = sName, x=DEgenes, rowNames = TRUE)
    saveWorkbook(workbook, 
                 file = paste0("DEgenes_", name, ".xlsx"),
                 overwrite = T)
  }
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<



# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Usage: fun.DEgenes_v2(seuratobject, assayname, samplename)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/DEanalysis_function.R
## version2 22/5/25
fun.DEgenes_v2 <- function(object, cluster, name){
  
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
    DEgenes$SYMBOL <- row.names(DEgenes)
    # 列の入れ替え
    DEgenes2　<- DEgenes[,c(6,1,2,3,4,5)]
    
    #シートへ保存：追加OK    
    sName = paste(stage[1],"vs", stage[j+1], sep="")
    addWorksheet(workbook, sheetName = sName)
    writeData(workbook, sheet = sName, x=DEgenes2, rowNames = FALSE)
    saveWorkbook(workbook, 
                 file = paste0("DEgenes_", cluster, "_", name, ".xlsx"),
                 overwrite = T)
  }
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## DEgenesを同定する = 条件で変動する遺伝子をFindMarkerで見つける
# Usage: function_DEanalysis.R(seuratobject, assayname, samplename)
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
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# Example : １度に全てのclusterについて解析を行いたい場合．
Idents(object) <- "seurat_clusters"
cluster.list　<- levels(object)

for(a_cluster in cluster.list){
  fun.DEgenes(object, a_cluster, "Mesenchyme1789Cells_SCT_Dim20Res0.2")
}

