## For GO
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(biomaRt)
library(GOSemSim)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(tidytree)


# #### Step3. Mouse databaseの準備　------------------------------------------
#マウスのデータベースはどこにある？
ensembl=useMart("ensembl")
listDatasets(ensembl)
## 107         mmusculus_gene_ensembl     Mouse genes (GRCm39)      GRCm39

db <- useMart(biomart = "ENSEMBL_MART_ENSEMBL") #使用するmart指定
Mm <- useDataset("mmusculus_gene_ensembl", mart = db)

# getBM()機能のfilters=引数、attributes=引数で指定可能なtermを調べるには、listFilters()機能、listAttributes()機能が有用である。

  ## attributeに何を入れる？
  listAttributes(Mm) ## 300くらい出てくる
  listAttributes(Mm)$name[grep(pattern = "symbol", listAttributes(Mm)$name)]
    # [1] "hgnc_symbol"       "mgi_symbol"        "uniprot_gn_symbol"
  listAttributes(Mm)$name[grep(pattern = "entrez", listAttributes(Mm)$name)]
    # [1] "entrezgene_trans_name"  "entrezgene_description" "entrezgene_accession"   "entrezgene_id"

  #  Keys <- rownames(Seurat_ojbjec)
  # 
  # all.genes.df <-
  #   getBM(attributes = c("mgi_symbol", "entrezgene_id"),      # ここで要素に入れたい項目を指定
  #         mart = Mm,
  #         filters = "mgi_symbol",
  #         values = Keys)


#> function --------------------------------------------
FUN.ClustProf.GOKGGRA.scRNAseq.mmu.v250118 <- function(
    df_DEG,
    Seurat_ojbject = NULL,
    ClusterofInterest,
    comparison2 = NULL,
    comparison3 = NULL,
    GeneColumn = "gene", 
    PadjColumn = "p_val_adj", 
    pathwiew_drawing = FALSE){
  
  # df_DEG <- df_Parent %>% filter(cluster == ClusterofInterest)
  # if(!exists(db, Mm, Keys, all.genes.df)) {
  print(" ------------------------------------------  ")
  print("Start ORA analysis in ")
  print(paste0(ClusterofInterest, comparison2, comparison3, sep = " "))
  print("Step1: Creating allGenes dataframe")
  
  # db <- useMart(biomart = "ENSEMBL_MART_ENSEMBL") #使用するmart指定
  # Mm <- useDataset("mmusculus_gene_ensembl", mart = db)
  
  Keys <- rownames(Seurat_ojbject)
  
  all.genes.df <-
    getBM(attributes = c("mgi_symbol", "entrezgene_id"),  #ここで要素に入れたい項目を指定
          mart = Mm,
          filters = "mgi_symbol",
          values = Keys)
  # } if exist
  
  # define file name
  filename <- ClusterofInterest
    if(is.null(comparison2) == FALSE){
        filename <- paste(filename, comparison2, sep="_")
        }
      if(is.null(comparison3) == FALSE){
          filename <- paste(filename, comparison3, sep="_")
          }
  new.dir <- paste0("CP_", filename)
  dir.create(new.dir)
  
  ## DEGのデータがある時のみ
  if(nrow(df_DEG) >0 ){
    
    #require(clusterProfiler)
    require(clusterProfiler)
    require(pathview)
    require(org.Mm.eg.db)
    ## 準備　--------------------------------------------------------
    ## Inputのdata frameにSymbol列を作成
    # names(df_DEG)[which(names(df_DEG) == GeneColumn)] <- "Symbol"
    # names(df_DEG)[which(names(df_DEG) == PadjColumn)] <- "p_val_adj"
    #　Inputの Symbolと referenceのmgi_symbolでくっつける
    
    ## coloumnのrenmae
    # df_DEG <- df_DEG %>% 
    #   dplyr::rename("gene" = GeneColumn)
    
    # data frameのmerge & entrezgene_idがSymbolとentrezgene_idのduplicationを除く
    df_DEG <- 
      df_DEG |>
        inner_join(all.genes.df, by = c("gene" = "mgi_symbol")) |>
        dplyr::filter(!is.na(entrezgene_id)) |>
        filter(!duplicated(entrezgene_id))
          print(head(df_DEG)) # output
      
    ## pvalue値にentrez idをつけたvector  
    genes_adjp <- df_DEG$p_val_adj
      names(genes_adjp) <- df_DEG$entrezgene_id 
    ## logFCへentrez idをつけた vector  
    Named_logFC_list <- df_DEG$avg_log2FC
      names(Named_logFC_list) <- df_DEG$entrezgene_id 
    
   print("Input genes are")
   print(genes_adjp)
      
   
    ##  2) GO  ############################# 
    res_GO <- data.frame()
    
    for(j in c("BP", "CC" ,"MF")){
      ego <- enrichGO(gene = names(Named_logFC_list), 
                      universe = all.genes.df$entrezgene_id,
                      OrgDb = org.Mm.eg.db,
                      ont = j, #vars
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.2,
                      minGSSize = 5,
                      readable = TRUE)
        # print(head(as.data.frame(ego)))
      
      if(nrow(as.data.frame(ego)) >= 1){  #結果が0なら飛ばす
        
        ego.simple <- simplify(ego)   # 似たような結果を除外する関数
        
        ego.simple.df <- ego.simple %>% 
          data.frame() %>% 
          dplyr::mutate(Ont = j)
        
        res_GO <- rbind(res_GO, ego.simple.df)
          # print(head(res_GO))
        
        GOtitle <- paste0("Enrich GO results of ", j," in ", filename)
        
        # 1 visualize the results
        p1 <- barplot(ego.simple, showCategory = 20) +
          ggtitle(GOtitle)
            plot(p1)
              ggsave(paste0(new.dir, "/1_GO_", j, "_dotplot", filename, ".pdf"),
                     width = 8, height = 12)
              
        #> ここからは結果が5つ以上の時にのみ表示する
        if(nrow(ego.simple)  > 5){
          # 2 visualize the dresults
          p2 <- clusterProfiler::cnetplot(ego.simple,
                                          showCategory = 5,
                                          node_label = "all",
                                          foldChange = Named_logFC_list) +
            ggtitle(GOtitle)
          
            plot(p2)
              ggsave(paste0(new.dir, "/1_GO_", j, "_centplot_", filename, ".png"),
                   width = 9 , height =12, dpi = 300)
          # 3 visualize the results
          ego_x2 <- enrichplot::pairwise_termsim(ego)
          
          p3 <- clusterProfiler::emapplot(ego_x2)  +
            ggtitle(GOtitle)
              plot(p3)
              ggsave(paste0(new.dir, "/1_GO_", j, "_emapplot_", filename, ".png"), 
                     width = 10, height =14)
           # if(j == "BP"){
              
          ## p4 <- treeplot(ego_x2)  +
              ##        ggtitle(GOtitle)
              ##       plot(p4)
            ##             ggsave(paste0(new.dir, "/1_GO_", j, "_Treeplot_", filename, ".png"),
              ##          width = 10, height =10)
            #    }
          # 4 visulalize GO DAG graph
          p5 <- goplot(ego.simple, showCategory = 10) +
            ggtitle(GOtitle)
            plot(p5)
              ggsave(paste0(new.dir, "/1_GO_", j, "_DAGgrph_", filename, ".png"),
                    width = 12, height =12)
              }# end resultが5以上の時
              
        print(paste0 ("GO:", j, " finished"))
        
      }# end if(nrow(as.data.frame(ego)) >=1)
    }# end loop j
    
    print("Finish all GO analysis")
 
    ## 2)  KEGG  ############################# 
    # http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
    print("Start KEGG analysis")
    
    df_kegg_null <- data.frame()
    df <- data.frame()
  
    ekg <- enrichKEGG(gene = names(Named_logFC_list), 
                      universe = all.genes.df$entrezgene_id,
                      organism = "mmu",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 5)
          # print("KEGG results are below: ")
          # print(head(as.data.frame(ekg)), 10)
    
    if(nrow(as.data.frame(ekg))>=1){
        ekg <- 
          ekg |> 
          setReadable(OrgDb = org.Mm.eg.db,
                      keyType = "ENTREZID")
      # mouseの記述を消す
      ekg@result$Description <- 
        gsub(pattern = " - Mus musculus (house mouse)",
             replacement = "", 
             ekg@result$Description, fixed = T)
      
      ekg@result <-
        ekg@result |>
        dplyr::filter(category != "Human Diseases") #human diseaseは除く
      
      # visualize the results
      pdf(paste0(new.dir, "/2_KEGG_bar_", filename, ".pdf"),
          width = 8, height = 12)
            barplot(ekg, showCategory = 10) + 
              ggtitle(paste0("Enrich KEGG results of ", filename))
                dev.off()
      # 2 visualize the results
      # p6 <- clusterProfiler::cnetplot(ekg,
      #                                 node_label = "all",
      #                                 categorySize = "pvalue",
      #                                 foldChange = Named_logFC_list)
      # plot(p6)
      # ggsave(paste0(new.dir, "/2_KEGG_centplot_", filename, ".png"),
      #        width = 10, height =12, dpi = 300)
      
      # # RESULT2
      res_KEGG <- as.data.frame(ekg)
      
      # if (nrow(res_KEGG)>=1){
      #   # pathwayによるpathwayデータを一括でかく
      ## logfcににentrez idをつけたvector  
      if(pathwiew_drawing == TRUE){
          genes_logfc <- df_DEG$avg_log2FC
          names(genes_logfc) <- df_DEG$entrezgene_id
          # print(head(genes_adjp))
          # print(head(genes_logfc))
        
          kegglist <- row.names(res_KEGG)
            for(i in 1:length(kegglist)){
              pathview(gene.data = genes_logfc, 
                        species = "mmu", 
                        pathway.id = kegglist[i],
                        limit = list(gene = 2, cpd = 1))
        }
      }
      
    }else{res_KEGG <- data.frame()}
    
    print("KEGG finished !")
    
    
    ## 3)  Reactome  ############################# 
    # https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
    print("start reactome")
    
    x <- enrichPathway(gene=names(Named_logFC_list),
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
    # print(as.data.frame(x))
    x_sum <- as.data.frame(x)
    
    if(nrow(x_sum>=1)){
      xb <- barplot(x) + labs(title = paste0("Reactome PA results of ", filename))
        plot(xb)
      ggsave(paste0(new.dir, "/3_RA_", filename, ".pdf"),
             width = 8, height = 12, dpi = 200)
    }else{x <- data.frame()}
    
    ## 4) dataのsave #########
    workbook <- createWorkbook()
      # sheetを作る
      addWorksheet(workbook, sheetName = "GO")
      addWorksheet(workbook, sheetName = "KEGG")
      addWorksheet(workbook, sheetName = "RA")
      # dataを書き込む
      writeData(workbook, sheet = 1, x = res_GO, rowNames = FALSE)
      writeData(workbook, sheet = 2, x = res_KEGG, rowNames = FALSE)
      writeData(workbook, sheet = 3, x = x_sum, rowNames = FALSE)
      #save
        saveWorkbook(workbook, 
                     file = paste0(new.dir, "/", filename, ".xlsx"), 
                      overwrite = TRUE)
    # リストにして返す
    saveRDS(list(NamedFoldChange = Named_logFC_list,
                 GO_res = ego, 
                 KEGG_res = ekg,
                 RA_res = x),
            paste0(new.dir, "/", filename, "_GOKEGGRAres.rds"))
    
  } # DEG_dfのif 終わり
  
  print("All ORA analysis finish !")
} # function終わり
## Functionここまで　-----------
