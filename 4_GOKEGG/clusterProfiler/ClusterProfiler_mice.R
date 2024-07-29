# Cluster profiler in zebrafish
library(clusterProfiler)
library(pathview) 
library(DOSE) 
library(AnnotationDbi)
library(org.Mm.eg.db)
# drawing
library(ggplot2)
library(ggrepel)
library(openxlsx)

# REFERENCES
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
# Q and A
# https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#faq
# (日本語) https://shiokoji11235.com/pathway-enrichment-analysis-in-r
## Cluster Profiler
######################### Description ###########################
# dataのinput
    # -> 有意な発現が見られる遺伝子リスト (FDR <0.01 など)：ENTREZIDで！
    # -> pvalueなどは必要ない
# https://shiokoji11235.com/go-analysis-for-transcriptome-analysis

# dataの準備：
data <- read.xlsx("/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp130_WH_scRNAseq_Analysis/CountCondition2(partialdata)/SeuratAnalysis/ExportFiles_Integration/220308_Analyzed_woRBCs/DEgenes_InterEpi_woRBC.61.15.xlsx", sheet=1)

# row.nameに1列目のsymbolを入れておく
row.names(data) <- data[,1]

#GENE SYMBOLやEnsembl ID をENTREZ IDへ変換し，結合する
id <- mapIds(
  org.Ms.eg.db, 
  keys=rownames(data),
  keytype="SYMBOL",
  column="ENTREZID")

# GENE SYMBOLをキーにして，dataframeへENTREZIDを追加する
data <- data.frame(data, ENTREZID=id) 

# ENTREZIDがNAのものを取り除く
  data <- data[!is.na(data$ENTREZID), ]
# ENTREZIDをrownameに切り変える
  row.names(data) <- data[,7]

  # adjp <0.01のリストを作成
  allgenes <- data[, 6] # vector型になる　is.vector(adjp)はTRUEが返ってくる．
  names(allgenes) <- data$ENTREZID # vectorへ名前を入れる : names関数
  adjp <- allgenes[allgenes<0.01] # adjusted pvalue (FDR) <0.01 だけ残す
  
  logfc <- data$avg_log2FC
  names(logfc) <- row.names(data)

  
## 1)  GO  ############################# 
ego <-
    enrichGO(
      gene = names(adjp),
      universe = row.names(data),
      OrgDb = org.Dr.eg.db,
      ont = "BP",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      minGSSize = 5,
      readable = TRUE)
head(as.data.frame(ego))

# 似たような結果を除外してくれるsimplifyという関数
ego.simple<-simplify(ego)
# 結果のsave:
write.xlsx(as.data.frame(ego.simple), "GO_ClusterProfiler.xlsx", rowNames = TRUE)
# visualize the results
clusterProfiler::dotplot(ego.simple)



## 2)  KEGG  ############################# 
# http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html

ekg <- 
  enrichKEGG(
    gene = names(adjp),
    universe = row.names(data),
    organism = "dre",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    minGSSize = 5)

head(as.data.frame(ekg))
browseKEGG(ekg,"dre04510")
browseKEGG(ekg,"dre04512")

# pathwayをみる
pathview(gene.data = logfc,
         species = "dre",
         pathway.id = "dre04510",
         limit = list(gene=2, cpd=1))


## 3)  GSEA  ##########################
## How to prepare gene list : https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist

gse <- gseGO(geneList = sort(allgenes, decreasing = TRUE),
             OrgDb = org.Dr.eg.db,
             keyType = "ENTREZID",
             ont = "BP",　　#"BP","CC","MF","ALL"
             minGSSize = 5,
             maxGSSize = 500,
             pAdjustMethod = "BH",
             pvalueCutoff = 0.05,
             eps = 0.0, 
             scoreType = "pos",
             verbose = FALSE)

# warnings
#2:  fgseaMultilevel(...) で: 
# For some pathways, in reality P-values are less than 1e-10. You can set the `eps` argument to zero for better estimation.
head((gse))
gseaplot(gse,geneSetID = "GO:0030198")




############ 　関数で一気に：有意なイベント出なくても空欄で出すようにした．　#########
########## GO, KEGG, GSEA function #################
fun.GOKGGSE <- function(data, Clustername, stage){
  # dataの作成：
  # row.nameに1列目のsymbolを入れておく
  row.names(data) <- data[,1]
  #GENE SYMBOLやEnsembl ID をENTREZ IDへ変換し，結合する
  id <- mapIds(org.Dr.eg.db, keys=rownames(data),
               keytype="SYMBOL", column="ENTREZID")
  # GENE SYMBOLをキーにして，dataframeへENTREZIDを追加する
  data <- data.frame(data, ENTREZID=id) 
  # ENTREZIDがNAのものを取り除く
  data <- data[!is.na(data$ENTREZID), ]
  # ENTREZIDをrownameに切り変える
  row.names(data) <- data[,7]
  
  # adjp, logFC のvectorを作成し，名前にentrezidを入れる
  allgenes <- data[, 6]
  names(allgenes) <- data$ENTREZID # vectorへ名前を入れる : names関数
  adjp <- allgenes[allgenes<0.01] # adjusted pvalue <0.01 を抽出する
  logfc <- data$avg_log2FC
  names(logfc) <- row.names(data)
  
  ## 1)  KEGG  ############################# 
  # http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
  ekg <- enrichKEGG(gene = names(adjp), 
                    universe = row.names(data),
                    organism = "dre",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 5)
  
  #readable trueが使えないので，ENTREZIDをsymbolへ変換する
  ekg <- setReadable(ekg, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
  # RESULT2
  res_KEGG <- as.data.frame(ekg)
  if (nrow(res_KEGG)>=1){
    # pathwayによるpathwayデータを一括でかく
    kegglist <- row.names(res_KEGG)
    
    for(i in 1:length(kegglist)){
      pathview(gene.data = logfc, 
               species = "dre", 
               pathway.id = kegglist[i],
               limit = list(gene=2, cpd=1))
    }
  }
  
  ##  2) GO  ############################# 
  ontology <- c("BP", "CC" ,"MF")
  res_GO <- data.frame()
  for(j in 1:length(ontology)){
    ego <- enrichGO(gene = names(adjp), 
                    universe = row.names(data),
                    OrgDb = org.Dr.eg.db,
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
      p1 <- clusterProfiler::dotplot(ego.simple)
      pdf(paste("GOdotplot",Clustername,stage, "_", ontology[j],".pdf", sep=""),
          width=8, height = 6)
      plot(p1)
      dev.off()
    }
  }
  
  ## 3)  GSEA  ##########################
  ## How to prepare gene list : https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist
  
  gse <- gseGO(geneList = sort(allgenes, decreasing = TRUE),
               OrgDb = org.Dr.eg.db,
               keyType = "ENTREZID",
               ont = "ALL",　　#"BP","CC","MF","ALL"
               minGSSize = 5,
               maxGSSize = 500,
               eps = 0, # 必要か？
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               scoreType = "pos",
               verbose = FALSE)
  # RESULT3
  
  gse <- setReadable(gse, OrgDb = org.Dr.eg.db, keyType = "ENTREZID")
  res_GSEA <- as.data.frame(gse)
  if (nrow(as.data.frame(gse))>=1){
    ID <- res_GSEA$ID
    Description <- res_GSEA$Description
    pdf(paste("GSEAplot", Clustername, stage, ".pdf", sep=""), width=8, height = 6)
    for(i in 1:length(ID)){
      g <- gseaplot(gse, geneSetID = ID[i], title = Description[i])
      plot(g)
    }
    dev.off()
  }
  
  
  ## 4) dataのsave ##########################
  filename <- paste("ClusterProfiler_", Clustername, "_", stage, sep="")
  workbook <- createWorkbook()
  # sheetを作る
  addWorksheet(workbook, sheetName = "GO")
  addWorksheet(workbook, sheetName = "KEGG")
  addWorksheet(workbook, sheetName = "GSEA")
  # dataを書き込む
  writeData(workbook, sheet = 1, x=res_GO, rowNames = FALSE)
  writeData(workbook, sheet = 2, x=res_KEGG, rowNames = FALSE)
  writeData(workbook, sheet = 3, x=res_GSEA, rowNames = FALSE)
  #save
  saveWorkbook(workbook, file = paste(filename,"xlsx", sep="."), overwrite = TRUE)
  
  # リストにしてかえす KEGG, GO, GESAの順
  l <- list(ekg, ego, gse)
  return(l)
}





#####  Install  #######################
# オブジェクト 'get_fun_from_pkg' がありません 
# 入らない時は，https://support.bioconductor.org/p/9139765/　を参照．rvcheckを入れるとううまくいくかも

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler", force=TRUE)
BiocManager::install("pathview", force=TRUE)





