########### Featureplot #######################################
# Geneを1つずつfeatureplotで書いて，まとめて保存する
# dependencyはpatchwork
fun.fp.batch <- function(object, genelist, samplename, row.num){
  ps <- list()
  for (gene in genelist){
    p <- FeaturePlot(object, gene)
    ps[[gene]] <- p
  }
  col.num <- ceiling(length(genelist)/row.num) #ceiling: x未満でない整数を返す
  wid <- 5.15*col.num
  hei <- 4.415*row.num
  wrap_plots(ps, nrow=row.num)
  
  ggsave(paste("FP_",samplename,".png",sep=""), dpi=300, width=wid, height=hei)
}
# ex)
  ECmarkers <- c("EGFP","kdrl","pecam1", "fli1a")
  fun.fp.batch(normalized.un, ECmarkers, "normalized.unwounded")



## FeaturePlot #######################################
# merge sampleをサンプルごとに描く
# https://divingintogeneticsandgenomics.rbind.io/post/customize-featureplot-in-seurat-for-multi-condition-comparisons-using-patchwork/
FeaturePlotSingle<- function(object, a_gene, metadata_column, ...){
  all_cells<- colnames(object)
  groups<- levels(object@meta.data[, metadata_column])
    
  # the minimal and maximal of the value to make the legend scale the same. 
  minimal <- min(object[['RNA']]@data[a_gene, ])
  maximal <- max(object[['RNA']]@data[a_gene, ])
  ps<- list()
  # メタデータのgroupごとにグラフをかく
  for (group in groups) {
    subset_indx<- object@meta.data[, metadata_column] == group
    subset_cells<- all_cells[subset_indx]
    p <- FeaturePlot(object, features = a_gene, cells= subset_cells, ...) +
         scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
         ggtitle(group) +
         theme(plot.title = element_text(size = 10, face = "bold"))
    ps[[group]] <- p
  }
  return(ps)
}  
    
p_list <- FeaturePlotSingle(pbmc, feature= "MS4A1", metadata_column = "samples", pt.size = 0.05, order =TRUE)




##########################################################
##### ２つのfeatureをmergeさせる : blendのoptionを追加するだけ
# https://staffblog.amelieff.jp/entry/2020/01/14/120000
fun.fp.blend.2genes <- function(object, two.genes, samplename, thsld){
  FeaturePlot(normalized.un,
              features = two.genes, 
              blend = TRUE,
              blend.threshold = thsld)
  pngname = paste("FPblend_",samplename,two.genes[1], two.genes[2],".png",sep="")
  ggsave(pngname, dpi=300, width=14, height=3.5)
}
rb <- c("mScarletI-C1-2A-epNTR","pdgfrb")
fun.fp.blend.2genes(normalized.un, rb, "normalized.unwounded", 0.2)




fun.cluster.condition <- function(obj, filename){
  DefaultAssay(obj) <- "integrated"
  
  obj <- RunPCA(obj, verbose = FALSE)
  
  pdfname <- paste("ClusterFinding-", filename, ".pdf", sep="")
  pdf(pdfname, width = 6, height = 6)
  
  for (i in 1:length(con_dim)){
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:con_dim[i])
    obj <- FindNeighbors(object = obj, dims = 1:con_dim[i])
    for (j in 1:length(con_res)){
      object <- FindClusters(object = obj, resolution = con_res[j])
      u <- UMAPPlot(object = obj , label=T) + 
        NoLegend() + 
        labs(title=paste("Dimention = ", con_dim[i], ",  Resolution = ", con_res[j]), sep="") + 
        theme(plot.title = element_text(hjust = 0.5))
      plot(u) 
    }
  }
  dev.off()
}


#### >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## DotPlot ########################






####################################################
############# dim plot #############################
fun.dimplot.int <- function(object, filename){
  
  pngprefix <- paste("Dimplot-", filename, sep="")

  #UMAP と sSNE plot
  DimPlot(object, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
    ggsave(paste(pngprefix, "_umap_cluster.png", sep=""), dpi=300, width=5, height=5)
  DimPlot(object, reduction = "tsne", label = TRUE, repel = TRUE) + NoLegend()
    ggsave(paste(pngprefix, "_tsne_cluster.png", sep=""), dpi=300, width=5, height=5)
  DimPlot(object, reduction = "umap", group.by = "Condition")
    ggsave(paste(pngprefix, "_umap_condition.png", sep=""), dpi=300, width=6, height=5)
  DimPlot(object, reduction = "tsne", group.by = "Condition")
    ggsave(paste(pngprefix, "_tsne_condition.png", sep=""), dpi=300, width=6, height=5)
  
  DimPlot(object, reduction = "umap", split.by = "Condition", label = TRUE)
    ggsave(paste(pngprefix, "_umap_each_dpi.png", sep=""), dpi=300, width=16, height=4)
  DimPlot(object, reduction = "tsne", split.by = "Condition", label = TRUE)
    ggsave(paste(pngprefix, "_tsne_each_dpi.png", sep=""), dpi=300, width=16, height=4)
}