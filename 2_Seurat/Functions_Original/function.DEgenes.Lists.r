## clusterごと，異なるステージのサンプルのDEGをそれぞれ比較し， excel にまとめていく

clusterlist <- c("Celltype1", "Celltype2", "Celltype3", "Celltype4", "Celltype5")
Stage <- c("Day0", "Day1", "Day3", "Day7")


fun.DEgenes <- function(object, cluster, name){
  
  Idents(object) <- "CellType_x_stage" ## すでにidentとしてあるものを指定
  #Cluster x stageでのobjectをそれぞれ作る
  group <- NULL #最初に作っておかないとエラーが出る
    for(i in 1:length(Stage)){
      group[i] <- paste(cluster, Stage[i], sep="_")
      }
  print(group)
  workbook <- createWorkbook()
  
  # Day0を元に 他のstageのDEGを解析する．
  for(j in 1:(length(group)-1)){
    
    DE <- FindMarkers(object,
                      ident.1 = group[j+1],
                      ident.2 = group[1],
                      logfc.threshold = 0.2,
                      min.pct = 0.20)
    
    DEgenes = DE %>% 
      filter(p_val_adj < 0.05)  %>% 
      mutate(Expression = case_when(
        avg_log2FC > 0 ~ "Up",
        avg_log2FC < 0 ~ "Down"
      ))
    
    #シートへ保存：追加OK    
    sName = paste(Stage[1], "vs", Stage[j+1], sep="")
    addWorksheet(workbook, sheetName = sName)
    writeData(workbook, sheet = sName, x = DEgenes, rowNames = TRUE)
    saveWorkbook(workbook, 
                 file = paste0("5_DEgenes_", name, "_",　cluster,  ".xlsx"),
                 overwrite = T)
  }
}


## Execution
for(a_cluster in cluster_list){

    fun.DEgenes(Object, a_cluster, "SampleName")

}
