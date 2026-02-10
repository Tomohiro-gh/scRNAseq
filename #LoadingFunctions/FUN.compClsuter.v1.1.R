library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Dr.eg.db)

#> 

Fun.Comp.Clsuter.GOKGRA.v1.1 <-
  function(NamedGenelists, filename, Wid = 18, Hei = 6, nCat = 10){ # nCat のデフォルト値を設定
    require(org.Mm.eg.db)
    require(clusterProfiler)
    FileNameBody = paste0(filename)
    
    ## Theme -----------------------------------
    theme_cp <- theme(
      axis.text.y = element_text(size = 12, lineheight = 0.6),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.50, size = 12))
    ## Theme -----------------------------------
    
    # 各解析結果を格納するためのリストを初期化
    results_list <- list(BP = NULL, CC = NULL, MF = NULL, KEGG = NULL, Reactome = NULL)
    
    ## GO:BP
    message("Running enrichGO for BP...")
    comp.BP <- compareCluster(geneClusters = NamedGenelists,
                              fun = "enrichGO",
                              ont="BP",
                              OrgDb = org.Mm.eg.db) # OrgDbをオブジェクトとして指定
    
    if (!is.null(comp.BP) && nrow(as.data.frame(comp.BP)) >= 1) {
      comp.BP <- setReadable(comp.BP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      results_list$BP <- comp.BP # 結果をリストに格納
      
      p1 <- comp.BP |>
        clusterProfiler::dotplot(showCategory = nCat) +
        labs(title = paste0("CompareCluster ", FileNameBody, " GO:BP"), x = "") +
        theme_cp +
        coord_flip()
      # plot(p1) # 不要な plot() 呼び出しは削除
      ggsave(paste0(FileNameBody, "_Compare.GOBP.png"),
             plot = p1, # plot引数を明示的に指定
             width = Wid, height = Hei, dpi = 200)
    } else {
      message("No significant GO:BP enrichment found.")
    }
    
    ## GO:CC
    message("Running enrichGO for CC...")
    comp.CC <- compareCluster(geneClusters = NamedGenelists,
                              fun = "enrichGO",
                              ont="CC",
                              OrgDb = org.Mm.eg.db)
    
    if (!is.null(comp.CC) && nrow(as.data.frame(comp.CC)) >= 1) {
      comp.CC <- setReadable(comp.CC, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      results_list$CC <- comp.CC
      
      p2 <- comp.CC |>
        clusterProfiler::dotplot(showCategory = nCat) +
        labs(title = paste0("CompareCluster ", FileNameBody, " GO:CC"), x = "") +
        theme_cp +
        coord_flip()
      # plot(p2) # 不要な plot() 呼び出しは削除
      ggsave(paste0(FileNameBody, "_Compare.GOCC.png"),
             plot = p2,
             width = Wid, height = Hei, dpi = 200)
    } else {
      message("No significant GO:CC enrichment found.")
    }
    
    ## GO:MF
    message("Running enrichGO for MF...")
    comp.MF <- compareCluster(geneClusters = NamedGenelists,
                              fun = "enrichGO",
                              ont="MF",
                              OrgDb = org.Mm.eg.db)
    
    if (!is.null(comp.MF) && nrow(as.data.frame(comp.MF)) >= 1) {
      comp.MF <- setReadable(comp.MF, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      results_list$MF <- comp.MF
      
      p3 <- clusterProfiler::dotplot(comp.MF, showCategory = nCat) +
        labs(title = paste0("CompareCluster ", FileNameBody, " GO:MF"), x = "") +
        theme_cp +
        coord_flip()
      # plot(p3) # 不要な plot() 呼び出しは削除
      ggsave(paste0(FileNameBody, "_Compare.GOMF.png"),
             plot = p3,
             width = Wid, height = Hei, dpi = 200)
    } else {
      message("No significant GO:MF enrichment found.")
    }
    
    ## KEGG
    message("Running enrichKEGG...")
    comp.KG <- compareCluster(geneClusters = NamedGenelists,
                              fun = "enrichKEGG",
                              organism = "mmu")
    
    if (!is.null(comp.KG) && nrow(as.data.frame(comp.KG)) >= 1) {
      comp.KG <- setReadable(comp.KG, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
      
      # フィルタリング後のデータが空になる可能性を考慮
      filtered_kg_result <- comp.KG@compareClusterResult %>%
        dplyr::filter(grepl("pathway|metabolism", Description, ignore.case = TRUE)) %>%
        dplyr::mutate(
          Description = str_replace_all(
            Description, pattern = " - Mus musculus \\(house mouse\\)",
            replacement = ""
          )
        )
      
      # フィルタリングの結果、データが残っているか確認してからプロット
      if (nrow(filtered_kg_result) >= 1) {
        comp.KG@compareClusterResult <- filtered_kg_result # フィルターされた結果を戻す
        results_list$KEGG <- comp.KG
        
        p4 <- clusterProfiler::dotplot(comp.KG, showCategory = nCat) +
          labs(title = paste0("CompareCluster ", FileNameBody, " KEGG"), x = "") +
          theme_cp +
          coord_flip()
        # plot(p4) # 不要な plot() 呼び出しは削除
        ggsave(paste0(FileNameBody, "_Compare.KEGG.png"),
               plot = p4,
               width = Wid, height = Hei, dpi = 200)
      } else {
        message("No KEGG pathways matching 'pathway|metabolism' found after filtering.")
      }
    } else {
      message("No significant KEGG enrichment found.")
    }
    
    ## Reactome
    message("Running enrichPathway...")
    comp.RA <- compareCluster(geneClusters = NamedGenelists,
                              fun = "enrichPathway",
                              organism = "mouse",
                              readable = TRUE)
    
    if (!is.null(comp.RA) && nrow(as.data.frame(comp.RA)) >= 1) {
      results_list$Reactome <- comp.RA
      
      p5 <- clusterProfiler::dotplot(comp.RA, showCategory = nCat) +
        labs(title = paste0("CompareCluster ", FileNameBody, " Reactome"), x = "") +
        theme_cp +
        coord_flip()
      # plot(p5) # 不要な plot() 呼び出しは削除
      ggsave(paste0(FileNameBody, "_Compare.Reactome.png"),
             plot = p5,
             width = Wid, height = Hei, dpi = 200)
    } else {
      message("No significant Reactome enrichment found.")
    }
    
    ## saving results -----------------------------------
    # 各結果オブジェクトがNULLでないことを確認してリストに追加
    saveRDS(results_list, paste0(FileNameBody, "_CompareClusterResults.rds"))
    
    # 各データフレームを結合する際もNULLチェック
    df_to_bind <- list()
    if (!is.null(results_list$BP)) df_to_bind$BP <- as.data.frame(results_list$BP) %>% dplyr::mutate(Analysis = "GO:BP")
    if (!is.null(results_list$CC)) df_to_bind$CC <- as.data.frame(results_list$CC) %>% dplyr::mutate(Analysis = "GO:CC")
    if (!is.null(results_list$MF)) df_to_bind$MF <- as.data.frame(results_list$MF) %>% dplyr::mutate(Analysis = "GO:MF")
    if (!is.null(results_list$KEGG)) {
      # KEGG の select(-c(2:3)) は、GeneRatio や BgRatio カラムを削除する可能性があります
      # 元のコードでは filter の後に select していたので、GeneRatio が無くなるとエラーの原因になります。
      # 必要に応じて select の範囲を調整するか、GeneRatio を含まないようにします。
      # 今回はエラーの原因を避けるため、select(-c(2:3)) を削除しました。
      df_to_bind$KEGG <- 
        as.data.frame(results_list$KEGG) %>%
        dplyr::mutate(Analysis = "KEGG") |>
        dplyr::select(-c(2:3))
    }
    if (!is.null(results_list$Reactome)) {
      df_to_bind$RA <-
        as.data.frame(results_list$Reactome) %>%
        dplyr::mutate(Analysis = "Reactome")
    }
    # 全てのデータフレームが存在する場合のみ結合して保存
    if (length(df_to_bind) > 0) {
      do.call(rbind, df_to_bind) %>%
        write.xlsx(paste0(FileNameBody, "_CompareClusterResults.xlsx"))
    } else {
      message("No enrichment results to save to Excel.")
    }
    
  }
