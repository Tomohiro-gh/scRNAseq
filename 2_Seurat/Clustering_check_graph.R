library(Seurat)
library(ggplot2)


## Unsuperviesed clustering 
## Clsuter annotationする前にcheckするfunction
fun.Clustering.Check.Graph <- 
  function(object, MyIdent1, MyIdent2, samplename){
    
    Idents(object) <- MyIdent1
    num <- nlevels(object)
    
    Idents(object) <- "seurat_clusters"
    
    name <- paste0("Unsuperviesed_", samplename)
    
    #1) vlnplot
    VlnPlot(object, 
            c("nFeature_RNA", "nCount_RNA","percent.mt"), 
            ncol=3, pt.size=0.1)
    ggsave(file = paste0(name, "_CountFeature_Vlnplot.png"),
           dpi = 300, width = 12, height = 4)
    
    #2) Featureplot
    FeaturePlot(object, 
                c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                ncol=3)
    ggsave(file = paste0(name, "_CountFeature_Featureplot.png"),
           dpi = 300, width = 20, height = 7)
    
    # 3) UMAPで確認
    # ident: seurat cluster
    UMAPPlot(object, 
             label=T)
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster.png"),
           dpi = 300, width = 7.7, height = 7) 
    
    UMAPPlot(object, 
             label=T) + 
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster_nolegend.png"),
           dpi = 300, width = 7, height = 7) 
    
    # ident: seurat cluster, split
    UMAPPlot(object, 
             split.by=MyIdent1, 
             label=T)
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster_split_",MyIdent1, ".png"),
           dpi = 300, width = 5*num, height = 5) 
    
    # ident: split ident
    UMAPPlot(object, 
             group.by=MyIdent1) 
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, ".png"),
           dpi = 300, width = 7, height = 7)
    # ident: split ident + split
    UMAPPlot(object, 
             group.by=MyIdent1,
             split.by=MyIdent1) + 
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1,"_split.png"),
           dpi = 300, width = 5*num, height = 5)
    
    #第二のidentを指定した場合は，下記を実行
    if(MyIdent2 !=""){
      #自分の見たいidentへ変更する
      UMAPPlot(object, 
               group.by=MyIdent2)
      ggsave(file = paste0(name, "_UMAP_by_,", MyIdent2,".png"),
             dpi = 300, width = 7.7, height = 7) 
      
      UMAPPlot(object, 
               group.by=MyIdent2,
               split.by=MyIdent2) + 
        NoLegend()
      ggsave(file = paste0(name, "_UMAP_by_", MyIdent2,"_split.png"),
             dpi = 300, width = 5*num, height = 5)
      
    }
  }


fun.Unspervised.Clustering.Check.Graph <- 
  function(
    SeuratObject, myCondition, myident, samplename){
  
  name <- paste0("Unsuperviesed_", samplename)
  
  #1) vlnplot
  VlnPlot(
    SeuratObject,
    c("nFeature_RNA", "nCount_RNA","percent.mt"),
    ncol=3, pt.size=0.1)
  ggsave(
    file = paste0(name, "_CountFeature_Vlnplot.png"),
    dpi = 300, width = 12, height = 4)
  
  #2) Featureplot
  FeaturePlot(
    SeuratObject,
    c("nFeature_RNA", "nCount_RNA","percent.mt"),
    ncol=3)
  ggsave(
    file = paste0(name, "_CountFeature_Featureplot.png"),
    dpi = 300, width = 20, height = 7)
  
  # 3) UMAPで確認
  Idents(SeuratObject) <- "seurat_clusters"
  # 3-1 Seurat clusterで
  UMAPPlot(
    SeuratObject, 
    label=T)
  ggsave(
    file = paste0(name, "_UMAP_byCluster.png"),
    dpi = 300, width = 7.5, height = 7) 
  
  # 3-2 指定したConditionでみる
  UMAPPlot(
    SeuratObject,
    group.by=myCondition) 
  ggsave(
    file = paste0(name, "_UMAP_byCondition_split.png"),
    dpi = 300, width = 7.5, height = 7)
  
  # 3-3 指定したConditionで splitして表示
  UMAPPlot(
    SeuratObject,
    split.by=myCondition,
    group.by=myCondition) +
    NoLegend()
  ggsave(
    file = paste0(name, "_UMAP_byCondition_split.png"),
    dpi = 300, width = 16, height = 4)
  
  # 3-4 さらに他にも見たいidentがあれば
  Idents(SeuratObject) <- myident
  UMAPPlot(
    SeuratObject, 
    group.by=myident) 
  ggsave(
    file = paste0(name, "_UMAP_by_",myident,".png"),
    dpi = 300, width = 7.5, height = 7)
  # 3-4 さらに indetでsplit
  UMAPPlot(
    SeuratObject, 
    group.by=myident,
    split.by=myident) +
    NoLegend()
  ggsave(
    file = paste0(name, "_UMAP_by_split_by_",myident,".png"),
    dpi = 300, width = 10, height = 6)
  
  
  }

fun.Unspervised.Clustering.Check.Graph(ECMC, "Condition", "State","ECMC_Dim3Res0.5")






## Cluster Annotation後の check graph
## Usage: 
## fun.Annotated.Clustering.Check.Graph(ECMC, "Annotation1", "State", "ECMC_dim3Res0.5")

fun.Annotated.Clustering.Check.Graph <- 
  function(
    SeuratObject, MyIdent, SplitCondition, samplename){
  
  name <- paste0("ClusterAnnotation_", samplename)
  
  Idents(SeuratObject) <- MyIdent
  # 1 Annotatedしたグラフで
  UMAPPlot(
    SeuratObject, label=T) + 
    labs(title="Cluster Annotation") + 
    theme(
      plot.title = 
        element_text(hjust=0.5, 
                     size=14, 
                     face="bold"),
      legend.position='none')
  ggsave(
    file = paste0(name, "_UMAP_label.png"),
    dpi = 300, width = 7, height = 7)
  
  # 2 
  UMAPPlot(
    SeuratObject) + 
    labs(title="Cluster Annotation") + 
    theme(
      plot.title = 
        element_text(hjust=0.5, 
                     size=14, 
                     face="bold"))
  ggsave(
    file = paste0(name, "_UMAP_legend.png"),
    dpi = 300, width = 7.5, height = 7) 
  
  # 3 Annotatedしたplotをconditionで分ける
  UMAPPlot(
    SeuratObject, 
    split.by=SplitCondition) + 
    labs(title="Cluster Annotation, split by condition") + 
    theme(
      plot.title = 
        element_text(hjust=0.5, 
                     size=14, 
                     face="bold"))
  ggsave(
    file = paste0(name, "_UMAP_split_by", SplitCondition,".png"),
    dpi = 300, width = 14, height = 7) 
  
}
