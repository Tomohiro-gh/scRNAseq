#### Seurat Graph section

fun.gencel2 <- function(data, filename, upper, lower){
  genes_per_cell <- Matrix::colSums(data>0)
  # count a gene only if it has non-zero reads mapped.
  
  filename <- as.character(filename)
  pngname <- paste("genes per cell (ordered)-", filename, ".png", sep="")
  
  png(pngname, width = 1000, height = 500)
  plot(sort(genes_per_cell), xlab='cell', log='y', main=paste(filename, 'genes per cell (ordered)', sep = "-"), cex = 0.6)
  abline(h=upper, col='red')  # upper threshold
  abline(h=lower, col='blue') # lower threshold
  dev.off()
}
#Execution example
# fun.gencel2(matrix.un, "1_unwounded", 4000, 70) # data, upper, lowerの順


## Dim Plot #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#### Dataの確認　: UMAPを使って
# Usage: fun.dimPlot.QC.umap(object, filename)
# /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/Seurat_Graph.R
fun.dimPlot.QC.umap <- function(object, filename){
  DimPlot(object, reduction = "umap", label = TRUE) + 
    NoLegend() +
    ggtitle(paste0("Dimplot ",filename)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Dimplot_", filename,".png"), 
         width=8, height=8, dpi=300)
  # sample条件ごとに色で分ける．1つのグラフ
  DimPlot(object, reduction = "umap", 
          group.by = "Condition", label = FALSE)
  ggsave(paste0("Dimplot_", filename,"_samplecolor.png"), 
         width=6, height=5, dpi=300)
  # グラフを分ける
  DimPlot(object, reduction = "umap", 
          split.by = "Condition", label = TRUE)
  ggsave(paste0("Dimplot_", filename,"_sampleSplit.png"), 
         width=16, height=4, dpi=300) 
}
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  
# Clustering前後の Dimplot, Featureplot, volcano plot for check
# Usage: fun.Clustering.Check.Graph(object, samplename,Condition,myident)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/Seurat_Graph.R
fun.Clustering.Check.Graph <- function(object, samplename, Condition,myident){
  # Dimplot
  # graph1
  DimPlot(object, reduction = "umap", label = TRUE) + 
    NoLegend() +
    ggtitle(samplename) +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    file = paste0(samplename, "_UMAP.png"), 
    dpi = 300, width = 7, height = 7)
  
  # vlnplot
  VlnPlot(
    object, 
    c("nFeature_RNA", "nCount_RNA","percent.mt"), 
    ncol=3, pt.size=0.1)
  ggsave(
    file = paste0(samplename, "_Vlnplot_CountFeature.png"), 
    dpi = 300, width = 15, height = 5)
  
  #FeaturePlot
  FeaturePlot(
    object,
    c("nFeature_RNA", "nCount_RNA","percent.mt"), 
    ncol=3)
  ggsave(
    file = paste0(samplename, "_Featureplot_CountFeature.png"),
    dpi = 300, width = 15, height = 5)
  
  # 3) UMAPで確認 # cluster情報で
  UMAPPlot(
    object, label=F)
  ggsave(
    file = paste0(samplename, "_UMAP_byCluster.png"),
    dpi = 300, width = 8, height = 7) 
  # Conditionごとに分ける
  UMAPPlot(
    object, 
    split.by=Condition, group.by=Condition) 
  ggsave(
    file = paste0(samplename, "_UMAP_byCondition_split.png"),
    dpi = 300, width = 16, height = 4)
  
  #自分の見たいidentへ変更する
  Idents(object) <- myident
  UMAPPlot(
    object, group.by=Condition) 
  ggsave(file = paste0(samplename, "_UMAP_by",Condition,".png"),
         dpi = 300, width = 8, height = 7)
}
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




################ Feature plot ################################
## ただplotするだけ．
fun.fp <- function(object, genes){
  fp <- FeaturePlot(object = object, 
                    features = genes, 
                    cols = c("grey", "blue"), 
                    reduction = "umap")
  plot(fp)
}
# Exmaples
# fun.fp(wh.int.38.15, "EGFP")


# 複数のgeneについて，Featureplotを描く　>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
########### Featureplot一気に表示して保存する
# Usage: fun.fp.batch(object, genelist, samplename, row.num)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/Seurat_Graph.R
fun.fp.batch <- function(object, genelist, samplename, row.num){
  require(patchwork)
  ps <- list()
  for (gene in genelist){
    
    p <- FeaturePlot(
      object, 
      gene,
      cols = c("#D5D5D3","#CB2314"))
    
    ps[[gene]] <- p
  }
  # geneリストの名前を全てのせる
  genenames <- paste(genelist, collapse = ",")
  col.num <- ceiling(length(genelist)/row.num) #ceiling: x未満でない整数を返す
  wid <- 5.15*col.num
  hei <- 4.415*row.num
  wrap_plots(ps, nrow=row.num)
  ggsave(paste0("FP_",samplename,"_",genenames,".png"), dpi=300, width=wid, height=hei)
}

# Example 1
#ECmarkers <- c("EGFP", "kdrl","pecam1","fli1a")
#fun.fp.batch(wh.int.sct.30.0.15, ECmarkers, "SCT30629Cells_Dim30Res0.15", 2)

# Exmaple2) リストになっている全てのオブジェクトについて一度に実行
#ECmarkers <- c("EGFP", "kdrl","pecam1","fli1a")
#for (i in 1:length(standard.4samples)){
#  a_object <- standard.4samples[[i]]
#  a_samplename <- names(standard.4samples)[i]
#  filname <- paste("standard.", a_samplename, ".ECmarkers", sep="")
#  fun.fp.batch(a_object, ECmarkers, filname, 2)
# }
# example 終了 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


## scoreingした後のグラフをかく．subtitleにどのgene
fun.fp.scored <- function(object, scorename, markers_vector ,samplename){
  require(wesanderson)
  FeaturePlot(
    object = object,
    features = scorename, 
    cols= wes_palette("Zissou1", 6, type="continuous")) +
    labs(subtitle = paste(markers_vector, collapse = ", ")) +
    theme(plot.subtitle = element_text(hjust=0.5))
  
  pngname = paste0("FP_scoring_", scorename,"_" ,samplename, ".png")
  ggsave(pngname, dpi=300, width=6, height=6)
  
}

  


################ Blendした図をかく ################################
fun.fp.blend.2genes <- function(object, two.genes, samplename, thsld){
  FeaturePlot(object,
              features = two.genes, 
              blend = TRUE,
              blend.threshold = thsld)
  pngname = paste("FPblend_", samplename, two.genes[1], two.genes[2],".png",sep="")
  ggsave(pngname, dpi=300, width=14, height=3.5)
}

## Example :
#rb <- c("mScarletI-C1-2A-epNTR","pdgfrb")
#for (i in 1:length(standard.4samples)){
#  a_object <- standard.4samples[[i]]
#  a_samplename <- names(standard.4samples)[i]
#  filname <- paste("standard.", a_samplename, sep="")
#  fun.fp.blend.2genes(a_object, rb, filname, 0.2)
# }
#ECs <- c("EGFP","kdrl")
#for (i in 1:length(standard.4samples)){
#  a_object <- standard.4samples[[i]]
#  a_samplename <- names(standard.4samples)[i]
#  filname <- paste("standard.", a_samplename, sep="")
#  fun.fp.blend.2genes(a_object, ECs, filname, 0.2)
#}



####################################################
################ Violin plot  ######################
fun.vp <- function(object, genes){
  vp <- VlnPlot(object, features = genes, split.by="Condition", pt.size = 0.1)
  plot(vp)
}
# fun.vp(wh.int.38.15, "cdh5")




## 複数の遺伝子を一気に書いて保存する
fun.vln.batch <- function(object, genelist, samplename){
  
  row.num <- length(genelist)
  ps <- list()
  for (gene in genelist){
    p <- VlnPlot(object, feature=gene, pt.size = 0.05) +
      NoLegend()
    ps[[gene]] <- p
  }
  col.num <- ceiling(length(genelist)/row.num) #ceiling: x未満でない整数を返す
  wid <- 10.3*col.num
  hei <- 2.2*row.num
  wrap_plots(ps, ncol=1)
  ggsave(paste0("Vlnplot_",samplename,".png"), dpi=300, width=wid, height=hei)
}

# Exmaple1) 
# ECmarkers <- c("EGFP", "kdrl","pecam1","fli1a")
# fun.vln.batch(wh.int.sct.30.0.15, ECmarkers, "SCT30629Cells_Dim30Res0.15")

# Exmaple2) リストになっている全てのオブジェクトについて一度に実行
# ECmarkers <- c("EGFP", "kdrl","pecam1","fli1a")
#for (i in 1:length(standard.4samples)){
#  a_object <- standard.4samples[[i]]
#  a_samplename <- names(standard.4samples)[i]
#  filname <- paste0("standard.", a_samplename, ".ECmarkers")
#  fun.vln.batch(a_object, ECmarkers, filname)
#}


##　identを指定して，特定のcluster における遺伝子発現
## nichenetrのtutorialより
fun.vln.selected <- function(
    SeuratObject, selectedident, genelist, splitCondition)
{
  row.num <- length(genelist)
  ps <- list()
  for (gene in genelist){
    p <- VlnPlot(
      SeuratObject %>% 
        subset(idents = selectedident),
      feature = gene, 
      split.by = splitCondition,
      pt.size = 0.05) +
      xlab("")#+NoLegend()
    ps[[gene]] <- p
  }
  
  col.num <- ceiling(length(genelist)/row.num) #ceiling: x未満でない整数を返す
  ## 横の長さを出す．
  yousosu <- SeuratObject@meta.data %>%
    select(splitCondition) %>%
    unique()
  wid <- 2*nrow(yousosu)
  hei <- 4*row.num
  wrap_plots(ps, ncol=1)
  ggsave(
    filename=
      paste0("Vlnplot_", splitCondition, "_in_", selectedident, "_in_",genelist[1],".png"),
    dpi=300, width=wid, height=hei)
}

## Example
#mygenes <- c("Robo4","Kdr")
#fun.vln.selected(Dim26Res0.5, "20", mygenes,"State")
#fun.vln.selected(Dim26Res0.5, "12", "Slit2","State")



fun.vln.selected.morethan.2groups <- function(
    SeuratObject, selectedident, genelist, GroupbyCondition)
{
  row.num <- length(genelist)
  ps <- list()
  for (gene in genelist){
    p <- VlnPlot(
      SeuratObject %>% 
        subset(idents = selectedident),
      feature = gene, 
      group.by = GroupbyCondition,
      pt.size = 0.05) +
      xlab("")#+NoLegend()
    ps[[gene]] <- p
  }
  
  col.num <- ceiling(length(genelist)/row.num) #ceiling: x未満でない整数を返す
  ## 横の長さを出す．
  yousosu <- SeuratObject@meta.data %>%
    select(GroupbyCondition) %>%
    unique()
  wid <- 2*nrow(yousosu)
  hei <- 4*row.num
  wrap_plots(ps, ncol=1)
  
  genenames = paste(genelist, collapse = ",")
  ggsave(
    filename=
      paste0("Vlnplot_", GroupbyCondition, "_in_", selectedident, "_in_", genenames, ".png"),
    dpi=300, width=wid, height=hei)
}


fun.vlp.split <- 
  function(seurat_object, features, group, split, name){
    VlnPlot(
      object = seurat_object,
      features = features,
      group.by=group,
      split.by=split,
      split.plot = TRUE,
      cols = c("#3B4992FF","#EE0000FF")) + 
      labs(x = "")
    
    ggsave(paste0("VlnPlot_",features, "_groupby_",group,"_split_by_",split, ".png"),
           width=8, height=4, dpi=300)
  }





################ Dot plot  ################################
fun.dp <- function(sub, markers, filename){
  pdfname <- paste0("DotplotEachCluster", filename, ".pdf")
  d <- DotPlot(sub, features=markers, cols=c("blue", "red","magenta", "cyan"), dot.scale=8, split.by="Condition") +
    RotatedAxis()
  pdf(pdfname, width=16, height=12)
  plot(d)
  dev.off()
}
# dotplot(wh.int.38.15, markers.dp, "Integration")


## DoHeatmap: scale dataを使用
# active identでgroupわけをする
fun.HM <- function(object, targetgenes)
  DoHeatmap(object = object,
            features = targetgenes,
            group.by = "ident") + 
  theme(axis.text.y = element_text(color = "black", size = 10))
# fun.HM(wh.int.woRBC.61.15, DEGs)

# DoHeatmapは，@ScaleData スロットを参照する
# defalutassayがRNAになっていると，scaledataのスロットはない．
# そのため，scale dataを作る
# https://satijalab.org/seurat/reference/doheatmap
# https://satijalab.org/seurat/archive/v3.1/merge_vignette.html
# https://github.com/satijalab/seurat/issues/2960
