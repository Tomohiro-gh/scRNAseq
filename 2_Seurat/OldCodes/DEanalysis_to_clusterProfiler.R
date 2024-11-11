library(Seurat)
library(ggplot2)
library(SeuratData)
library(patchwork)
library(Matrix)
library(gdata)
library(reshape)
library(reshape2)
library(openxlsx)
library(tidyverse)
library(clusterProfiler)
library(pathview) 
library(DOSE) #for GSEA
library(AnnotationDbi)
library(org.Dr.eg.db)



#### Seurat find markers functionなどで出した degを元に，

#### clusterProfiler (GO, KEGG, GSEA)を実施．

## Preparation :
## Allgene（今回対象とする全遺伝子）のobjectを作成する
allgene <- rownames(SCT22375[["RNA"]])
length(allgene)  #23158
allgenes.df <- data.frame(allgene, row.names = allgene)
id <- mapIds(org.Dr.eg.db, keys=rownames(allgenes.df),
             keytype="SYMBOL", column="ENTREZID")
id2 <- mapIds(org.Dr.eg.db, keys=rownames(allgenes.df),
              keytype="SYMBOL", column="ENSEMBL")
genes.df <- data.frame(ENTREZID=id, ENSEMBL=id2)
genes.df$SYMBOL <- row.names(genes.df)

# Merge result with ids      
allgenes.df <- merge(allgenes.df, genes.df, by="row.names", all=T)
allgenes.df <- allgenes.df[, -(1:2)]




### 22/06/01 modified -> 06/02/22 modified
########## GO, KEGG, GSEA function #################
# Usage: fun.clusterProfiler.GOKGGSE(data, Clustername, stage, returnObject)
# "/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/DEanalysis_to_clusterProfiler.R"

fun.clusterProfiler.GOKGGSE <- function(data, Clustername, stage, returnObject){
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
  data.sig <- data[data$p_val_adj<0.01, ]
  
  genes_adjp <- data.sig[, 6]
  names(genes_adjp) <- data.sig$ENTREZID # vectorへ名前を入れる : names関数
  #adjp <- allgenes[allgenes<0.01] # adjusted pvalue <0.01 を抽出する
  genes_logfc <- data.sig$avg_log2FC
  names(genes_logfc) <- row.names(data.sig)
  #確認
  print(head(genes_adjp), 10)
  
  ## 1)  KEGG  ############################# 
  # http://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
  ekg <- enrichKEGG(gene = names(genes_adjp), 
                    universe = allgenes.df$ENTREZID,
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
      pathview(gene.data = genes_logfc, 
               species = "dre", 
               pathway.id = kegglist[i],
               limit = list(gene=2, cpd=1))
    }
  }
  
  ##  2) GO  ############################# 
  ontology <- c("BP", "CC" ,"MF")
  res_GO <- data.frame()
  for(j in 1:length(ontology)){
    ego <- enrichGO(gene = names(genes_adjp), 
                    universe = allgenes.df$ENTREZID,
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
      pdf(paste0("GOdotplot_",Clustername,"_", stage, "_", ontology[j],".pdf"),
          width=8, height = 6)
      plot(p1)
      dev.off()
    }
  }
  
  ## 3)  GSEA  ##########################
  ## How to prepare gene list : https://yulab-smu.top/biomedical-knowledge-mining-book/faq.html#genelist
  
  gse <- gseGO(geneList = sort(genes_adjp, decreasing = TRUE),
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
    pdf(paste0("GSEAplot", Clustername, stage, ".pdf"), width=8, height = 6)
    for(i in 1:length(ID)){
      g <- gseaplot(gse, geneSetID = ID[i], title = Description[i])
      plot(g)
    }
    dev.off()
  }
  
  ## 4) dataのsave ##########################
  filename <- paste0("ClusterProfiler_", Clustername, "_", stage)
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
  
  # リストにしてかえすか，引数で yesを選択すれば 帰ってくる
  if(returnObject == "YES"){
    l <- list(ekg, ego, gse)
    return(l)
  }
  
}

# listとして返す場合
cp_C93 <- fun.clusterProfiler.GOKGGSE(C93, cluster.list[10], shtname[3], "YES")  
# listとして返さない場合
fun.clusterProfiler.GOKGGSE(C93, cluster.list[10], shtname[3], "NO")  
