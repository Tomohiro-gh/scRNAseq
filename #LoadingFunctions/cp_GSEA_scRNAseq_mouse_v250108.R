## For GO
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(biomaRt)
library(GOSemSim)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)


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

# Keys <- rownames(Seurat_ojbjec)
# 
# all.genes.df <-
#   getBM(attributes = c("mgi_symbol", "entrezgene_id"),  #ここで要素に入れたい項目を指定
#         mart = Mm,
#         filters = "mgi_symbol",
#         values = Keys)
  
  ## GO KEGG analysis Function ------------------------------------------------------------------------------------
FUN.cp.GOKGGRA.scRNAseq.mmu.v1 <- function(
    Seurat_ojbject,
    df_DEG, 
    ClusterofInterest,
    comparison2 = NULL,
    comparison3 = NULL,
    GeneColumn = "gene", 
    PadjColumn = "p_val_adj", 
    pathwiew_drawing = FALSE){
  
  # df_DEG <- df_Parent %>% filter(cluster == ClusterofInterest)
# if(!exists(db, Mm, Keys, all.genes.df)) {
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
    if(!is.null(comparison2) == FALSE){
      filename <- paste(filename, comparison2, sep="_")
    }
    if(is.null(comparison3) == FALSE){
      filename <- paste(filename, comparison3, sep="_")
    }
  new.dir <- paste0("CP_", filename)
    dir.create(new.dir)
  print(paste0("Filename is ", filename))
  
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
  
    df_DEG <- inner_join(df_DEG, all.genes.df, by = c("gene" = "mgi_symbol"))
      dim(df_DEG) # [1] 769   8
      print(head(df_DEG))
    
    # entrezgene_idがSymbolとentrezgene_idのduplicationを除く
    df_DEG <- df_DEG %>%  
      filter(!is.na(entrezgene_id)) %>% 
      filter(!duplicated(entrezgene_id))
    ## pvalue値にentrez idをつけたvector  
    genes_adjp <- df_DEG$p_val_adj
      names(genes_adjp) <- df_DEG$entrezgene_id # vectorへ名前を入れる : names関数
    ## logFCへentrez idをつけたvector  
    Named_logFC_list <- df_DEG$avg_log2FC
      names(Named_logFC_list) <- df_DEG$entrezgene_id # vectorへ名前を入れる : names関
  
    
 
  # 2)  GSEA  ##########################
  # How to prepare gene list : https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist
  # Fullgenesetが必要．
  FM_Allgenes <- FindMarkers(
    Seurat_ojbject,
    # group.by = "EC_int_cca_sub",
    ident.1 = ClusterofInterest, #active identである必要がある
    # ident.2 = ,
    min.cells.group = 1,
    min.cells.feature = 1,
    min.pct = 0,
    logfc.threshold = 0,
    only.pos = FALSE) %>%
    mutate(Symbol = rownames(.)) %>%
    inner_join(., all.genes.df, by = c("Symbol" = "mgi_symbol")) %>%
    #　Inputの Symbolと referenceのmgi_symbolでくっつける    #> Gene idの変換
    filter(!is.na(entrezgene_id)) %>%
    filter(!duplicated(entrezgene_id))
        # entrezgene_idがSymbolとentrezgene_idのduplicationを除く

  print("Start GSEA analysis")

  allgenes_logFC <- FM_Allgenes$avg_log2FC
  names(allgenes_logFC) <- rownames(FM_Allgenes)

  gse <- gseGO(geneList = sort(allgenes_logFC, decreasing = TRUE),
                OrgDb = org.Mm.eg.db,
                keyType = "ENTREZID",
                ont = "BP",　　#"BP","CC","MF","ALL"
                minGSSize = 5,
                maxGSSize = 500,
                eps = 0, # 必要か？
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                scoreType = "pos",
                verbose = FALSE)

  ## GSEAの結果をdotplotで
  nTop = 20
  d <- clusterProfiler::dotplot(gse, showCategory = nTop) +
    ggtitle(paste0("GSEA result Top", nTop," in ", ClusterofInterest))
  plot(d)
  ggsave(paste0(new.dir, "/2_dotplot_GSEA_", ClusterofInterest, ".pdf"),
         width = 8, height = 12, dpi = 300)

  if (nrow(as.data.frame(gse)) >= 1){

  ## data frameへ変換
  res_GSEA <- setReadable(gse,
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID") %>%
    as.data.frame()

  res_GSEA_selected <- res_GSEA %>%
    #  filter(qvalue == "BP") %>%
    filter(qvalue < 1e-5) %>%
    filter(setSize >= 10 & setSize <= 200)

    # ID <- res_GSEA$ID
    # Description <- res_GSEA$Description
    if (nrow(res_GSEA_selected) > 20){
      number = 5
    }else{number = nrow(res_GSEA_selected)}

    for(j in 1:number){
      g <- gseaplot(gse,
                    by = "runningScore",
                    geneSetID = res_GSEA_selected$ID[j],
                    title = res_GSEA_selected$Description[j])
        plot(g)
          # save
          ggsave(paste0(new.dir, "/2_GSEAplot_", ClusterofInterest, "_", res_GSEA_selected$Description[j],".png"),
                  width = 8, height = 6, dpi = 300)
    }
  }else{res_GSEA <- data.frame()}

  print("GSEA finished !")

  ## 3)  KEGG  ############################# 
    # http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
    ekg <- enrichKEGG(gene = names(Named_logFC_list), 
                      universe = all.genes.df$entrezgene_id,
                      organism = "mmu",
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "BH",
                      minGSSize = 5)
    print("KEGG results are below: ")
    print(as.data.frame(ekg))
    
    if(nrow(as.data.frame(ekg))){
      ekg <- setReadable(ekg, 
                         OrgDb = org.Mm.eg.db,
                         keyType = "ENTREZID")
      # visualize the results
      p5 <- barplot(ekg, showCategory = 10)
          pdf(paste0(new.dir, "/3_KEGG_bar_", filename, ".pdf"),
              width = 8, height = 12)
            plot(p5)
              dev.off()
      # 2 visualize the results
      p6 <- clusterProfiler::cnetplot(ekg,
                                      node_label = "all",
                                      categorySize = "pvalue",
                                      foldChange = Named_logFC_list)
        plot(p6)
          ggsave(paste0(new.dir, "/3_KEGG_centplot_", filename, ".png"),
           width = 10, height =12, dpi = 300)
          
    # # RESULT2
    res_KEGG <- as.data.frame(ekg)
    
    }else{res_KEGG <- data.frame()}
   
   print("KEGG finished !")
    
    ## 4)  Reactome  ############################# 
    # https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
    print("start reactome GSE")
    x <- enrichPathway(gene=names(Named_logFC_list),
                       organism = "mouse",
                       pvalueCutoff = 0.05,
                       readable = TRUE)
    print(as.data.frame(x))
    x_sum <- as.data.frame(x)
    
    if(nrow(x_sum>=1)){
     xb <- barplot(x) + labs(title = "Reactome PA")
     plot(xb)
        ggsave(paste0(new.dir, "/3_RA_", filename, ".png"),
               width = 8, height = 10, dpi = 300)
    }else{x <- data.frame()}

    # y <- gsePathway(geneList = sort(allgenes_logFC, decreasing = TRUE),
    #                 organism = "mouse",
    #                 minGSSize = 5,
    #                 verbose = FALSE)
    # y_sum <- as.data.frame(y)
    # if(!is.null(y)){
    #   clusterProfiler::dotplot(y) + labs(title = "gsePathway (RA)")
    #     ggsave(paste0(new.dir, "/4_RAgesa_", ClusterofInterest, ".png"),
    #            width = 10, height =12, dpi = 300)
    # }else{y <- data.frame()}
    # 
    ## 5) dataのsave #########

    workbook <- createWorkbook()
    # sheetを作る
    addWorksheet(workbook, sheetName = "GO")
    #addWorksheet(workbook, sheetName = "GOGSEA")
    addWorksheet(workbook, sheetName = "KEGG")
    addWorksheet(workbook, sheetName = "RA")
    # addWorksheet(workbook, sheetName = "RAGSEA")
    # dataを書き込む
    writeData(workbook, sheet = 1, x = res_GO, rowNames = FALSE)
    #writeData(workbook, sheet = 2, x = res_GSEA, rowNames = FALSE)
    writeData(workbook, sheet = 2, x = res_KEGG, rowNames = FALSE)
    writeData(workbook, sheet = 3, x = x_sum, rowNames = FALSE)
    #writeData(workbook, sheet = 5, x = y_sum, rowNames = FALSE)
    #save
    saveWorkbook(workbook, 
                 file = paste0(new.dir, "/", filename, ".xlsx"), 
                 overwrite = TRUE)
    # リストにしてかえすか，引数で yesを選択すれば 帰ってくる
    # if(returnObject == "YES"){
    saveRDS(list(NamedFoldChange = Named_logFC_list,
                 GO_res = ego, 
                 KEGG_res = ekg,
                 RA_res = x),
            paste0(new.dir, "/", filename,"_GOKEGGRAres.rds"))
     
    } # DEG_dfのif 終わり
  } # function終わり
## Functionここまで　-----------
