## 
## DEgenetable:
## Allgenes: ENTREZIDで準備 rownames(allgens)などでもいい

fun.clusterProfiler.GO.KEGG.Human <-
  function(DEgeneTable, AllGenes, saveDir, samplename){
    require(clusterProfiler)
    require(org.Hs.eg.db)
    
  wd = getwd()
  
  dir.create(saveDir)
  setwd(paste0("./", saveDir))
    
  id <- mapIds(
    org.Hs.eg.db, 
    keys=rownames(DEgeneTable),
    keytype="SYMBOL",
    column="ENTREZID")
  # GENE SYMBOLをキーにして，dataframeへENTREZIDを追加する
  DEgeneTable <- data.frame(DEgeneTable, ENTREZID=id) 
  # ENTREZIDがNAのものを取り除く
  DEgeneTable <- DEgeneTable[!is.na(DEgeneTable$ENTREZID), ]
  # ENTREZIDをrownameに切り変える
  row.names(DEgeneTable) <- DEgeneTable$ENTREZID
  
  ## GO anaysis  
  ontology <- c("BP", "CC" ,"MF")
  res_GO <- data.frame()
    
# BP, CC, MFの３つを行う
for(a_ont in ontology){
  ego <-
    enrichGO(
      gene = row.names(DEgeneTable),
      universe = AllGenes,
      OrgDb = org.Hs.eg.db,
      ont = a_ont, #vars
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      minGSSize = 5,
      readable = TRUE)
      
  # 似たような結果を除外してくれるsimplifyという関数
  ego.simple <- simplify(ego)
  res <- as.data.frame(ego.simple)
      
  if (nrow(as.data.frame(res))>=1){  
      res$ont = a_ont
      res_GO <- rbind(res_GO, res)
        
  # visualize the results
  p <- clusterProfiler::dotplot(
    ego.simple,
    showCategory=20)
  p
  ggsave(
    paste0(samplename, "_GO_", a_ont, "_top20.png"),
    width=7, height = 10, dpi=300)
    }
  }
  # 結果のsave:
  write.xlsx(
    res_GO,
    paste0(samplename, "_GO_clusterProfiler.xlsx"),
    rowNames = FALSE)
    
## 2)  KEGG  ############################# 
ekg <- 
  enrichKEGG(
    gene = row.names(DEgeneTable),
    universe = AllGenes,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5)
  # data save
  write.xlsx(
    as.data.frame(ekg), 
    paste0(samplename, "_KEGG_ClusterProfiler.xlsx"), 
    rowNames = FALSE)
  # graph save
  p <- clusterProfiler::dotplot(ekg)
  p
  ggsave(
    paste0(samplename, "_KEGG_clusterProfiler_top20.png"),
    width=9, height=9, dpi = 300)
    
　# pathway 全部抽出
  DE_pathways <-
    as.data.frame(ekg) %>% 
    rownames()
    
  logfc <- DEgeneTable$avg_log2FC
  names(logfc) <- row.names(DEgeneTable)
    
  pathviewDir = paste0("pathview")
  dir.create(pathviewDir)
  setwd(pathviewDir)
  
  for(a_pathway in DE_pathways){
    pathview(
      gene.data = logfc,
      species = "hsa",
      pathway.id = a_pathway,
      limit = list(gene=2, cpd=1))}
  
  # 元のwdへ戻る
  setwd(wd)
}