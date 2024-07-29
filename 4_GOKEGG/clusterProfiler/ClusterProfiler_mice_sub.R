fun.clusterProfiler.GOKG.mmu <- function(df_DEG, Clustername, stage, UporDown){
  new.dir <- paste0("6_GOKEGG_", Clustername, "_3monthvs", stage, "_", UporDown)
  dir.create(new.dir)
  
  #require(clusterProfiler)
  require(clusterProfiler)
  require(pathview)
  ## 準備　--------------------------------------------------------
  ## Inputのdata frameにSymbol列を作成
  df_DEG$Symbol <- rownames(df_DEG)
  #　Inputの Symbolと referenceのmgi_symbolでくっつける
  df_DEG <- inner_join(df_DEG, all.genes.df, by= c("Symbol" = "mgi_symbol"))
  dim(df_DEG) # [1] 769   8
  print(head(df_DEG))
  
  # entrezgene_idがSymbolとentrezgene_idのduplicationを除く
  df_DEG <- df_DEG %>%  
    filter(!is.na(entrezgene_id)) %>% 
    filter(!duplicated(entrezgene_id))
  dim(df_DEG) # [1] 762   8
  
  ## pvalue値にentrez idをつけたvector  
  genes_adjp <- df_DEG$p_val_adj
  names(genes_adjp) <- df_DEG$entrezgene_id # vectorへ名前を入れる : names関数
  ## logfcににentrez idをつけたvector  
  genes_logfc <- df_DEG$avg_log2FC
  names(genes_logfc) <- df_DEG$entrezgene_id
  # print(head(genes_adjp))
  # print(head(genes_logfc))
  
  
  ## 1)  KEGG  ############################# 
  # http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
  ekg <- enrichKEGG(gene = names(genes_adjp), 
                    universe = all.genes.df$entrezgene_id,
                    organism = "mmu",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 5)
  #readable trueが使えないので，ENTREZIDをsymbolへ変換する
  ekg <- setReadable(ekg, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  
  # visualize the results
  p <- clusterProfiler::dotplot(ekg, showCategory = 20)
  pdf(paste0(new.dir, "/6_KGdotplot_", Clustername, "_", stage, "_", UporDown, ".pdf"),
      width = 8, height = 12)
  plot(p)
  dev.off()
  # # RESULT2
  res_KEGG <- as.data.frame(ekg)
  # if (nrow(res_KEGG)>=1){
  #   # pathwayによるpathwayデータを一括でかく
  #   kegglist <- row.names(res_KEGG)
  #   
  #   for(i in 1:length(kegglist)){
  #     pathview(gene.data = genes_logfc, 
  #              species = "mmu", 
  #              pathway.id = kegglist[i],
  #              limit = list(gene=2, cpd=1))
  #   }
  # }
  # 
  ##  2) GO  ############################# 
  ontology <- c("BP", "CC" ,"MF")
  res_GO <- data.frame()
  for(j in 1:length(ontology)){
    ego <- enrichGO(gene = names(genes_adjp), 
                    universe = all.genes.df$entrezgene_id,
                    OrgDb = org.Mm.eg.db,
                    ont = ontology[j], #vars
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2,
                    minGSSize = 5,
                    readable = TRUE)
    # 似たような結果を除外してくれるsimplifyという関数
    if (nrow(as.data.frame(ego))>=1){
      ego.simple <- simplify(ego)
      res_GO <- as.data.frame(ego.simple)
      res_GO <- rbind(res_GO, as.data.frame(ego.simple))
      
      # visualize the results
      p1 <- clusterProfiler::dotplot(ego.simple, showCategory = 20)
      pdf(paste0(new.dir, "/6_GOdotplot_",Clustername,"_", stage, "_",  UporDown, "_",  ontology[j],".pdf"),
          width = 8, height = 12)
      plot(p1)
      dev.off()
    }
  }
  
  ## 3)  GSEA  ##########################
  ## How to prepare gene list : https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist
  
  # gse <- gseGO(geneList = sort(genes_adjp, decreasing = TRUE),
  #              OrgDb = org.Mm.eg.db,
  #              keyType = "ENTREZID",
  #              ont = "ALL",　　#"BP","CC","MF","ALL"
  #              minGSSize = 5,
  #              maxGSSize = 500,
  #              eps = 0, # 必要か？
  #              pAdjustMethod = "BH",
  #              pvalueCutoff = 0.05,
  #              scoreType = "pos",
  #              verbose = FALSE)
  # # 
  # gse <- setReadable(gse, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
  # 
  # res_GSEA <- as.data.frame(gse)
  # 
  # if (nrow(res_GSEA)>=1){
  #   ID <- res_GSEA$ID
  #   Description <- res_GSEA$Description
  #   pdf(paste0("6_GSEAplot", Clustername, stage, ".pdf"), width=8, height = 6)
  #     for(i in 1:length(ID)){
  #       g <- gseaplot(gse, geneSetID = ID[i], title = Description[i])
  #         plot(g)
  #   }
  #   dev.off()
  # }
  # 
  ## 4) dataのsave ##########################
  filename <- paste0("6_ClusterProfiler_", Clustername, "_", stage)
  workbook <- createWorkbook()
  # sheetを作る
  addWorksheet(workbook, sheetName = "GO")
  addWorksheet(workbook, sheetName = "KEGG")
  # addWorksheet(workbook, sheetName = "GSEA")
  # dataを書き込む
  writeData(workbook, sheet = 1, x=res_GO, rowNames = FALSE)
  writeData(workbook, sheet = 2, x=res_KEGG, rowNames = FALSE)
  # writeData(workbook, sheet = 3, x=res_GSEA, rowNames = FALSE)
  #save
  saveWorkbook(workbook, file = paste0(new.dir, "/", filename, "_", UporDown ,".xlsx"), overwrite = TRUE)
  
  # リストにしてかえすか，引数で yesを選択すれば 帰ってくる
  # if(returnObject == "YES"){
  l <- list(ekg, ego)
  return(l)
  #} 
}
## ------------------------------------------------------------------------------------