library(Seurat)
library(ggplot2)


## Unsuperviesed clustering 
## Clsuter annotationする前にcheckするfunction
# fun.Unspervised.Clustering.Check.Graph(object, MyIdent1, MyIdent2, samplename)
fun.Unspervised.Clustering.Check.Graph <- 
  function(object, MyIdent1, MyIdent2, samplename){
    
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
    Idents(object) <- "seurat_clusters"
    
    UMAPPlot(object, 
             label=T) + 
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster.png"),
           dpi = 300, width = 7, height = 7) 
    
    UMAPPlot(object, 
             label=F) 
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster_legend.png"),
           dpi = 300, width = 7.7, height = 7) 
    
  
    
    Idents(object) <- MyIdent1
    # Idents(object) <- "seurat_clusters"
    num = Idents(object) %>% unique %>% length()
    
    UMAPPlot(object, 
             label=F)
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, ".png"),
           dpi = 300, width = 7, height = 7) 
    # ident: group by ident
    UMAPPlot(object, 
             group.by=MyIdent1,
             label=F) 
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, ".png"),
           dpi = 300, width = 7, height = 7)
    
    UMAPPlot(object, 
             group.by=MyIdent1, 
             split.by=MyIdent1,
             label=F) + 
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, "_split.png"),
           dpi = 300, width = 5*num, height = 7)
    # clusterring を groupでsplit
    Idents(object) <- "seurat_clusters"
    UMAPPlot(object, 
             label=F,
             split.by = MyIdent1)  + 
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster_split_by_", MyIdent1, ".png"),
           dpi = 300, width = 5*num, height = 7) 
    
    
    #第二のidentを指定した場合は，下記を実行
    if(MyIdent2 !=""){
      Idents(object) <- MyIdent2
      num = Idents(object) %>% unique %>% length()
      
      #自分の見たいidentへ変更する
      UMAPPlot(object, 
               split.by=MyIdent2)
      ggsave(file = paste0(name, "_UMAP_by_split_by_", MyIdent2,".png"),
             dpi = 300, width = 5*num, height = 5) 
      
      UMAPPlot(object, 
               group.by=MyIdent2)
      ggsave(file = paste0(name, "_UMAP_by_", MyIdent2,".png"),
             dpi = 300, width = 7.7, height = 7) 
      
      UMAPPlot(object, 
               group.by=MyIdent2,
               split.by=MyIdent2) + 
        NoLegend()
      ggsave(file = paste0(name, "_UMAP_by_", MyIdent2,"_split.png"),
             dpi = 300, width = 5*num, height = 5)
      
    }
  }


## Unsuperviesed clustering 
## Clsuter annotationする前にcheckするfunction
# fun.Unspervised.Clustering.Check.Graph.v2(object, MyIdent1, MyIdent2, samplename)
fun.Unspervised.Clustering.Check.Graph.v2 <- function(
    object, MyIdent1, MyIdent2, samplename){
  
  name <- paste0("Unsuperviesed_", samplename)
  
  #1) vlnplot
  VlnPlot(object, c("nFeature_RNA", "nCount_RNA", 
                    "percent.mt", "percent.ribo",
                    "S.Score", "G2M.Score"), 
          ncol = 2, pt.size = 0.1)
  ggsave(file = paste0(name, "_CountFeature_Vlnplot.png"), 
         dpi = 300, width = 8, height = 16)
  
  #2) FeaturePlot
  FeaturePlot(object, c("nFeature_RNA", "nCount_RNA", 
                        "percent.mt", "percent.ribo",
                        "S.Score", "G2M.Score"), ncol = 2)
  ggsave(file = paste0(name, "_CountFeature_Featureplot.png"), 
         dpi = 300, width = 8, height = 12)
  
  # 3) UMAPで確認
  # ident: seurat cluster
  # Idents(object) <- "seurat_clusters"
  
  UMAPPlot(object, label = T) + 
    NoLegend() + 
    labs(title="Unsupervised Clustering") + 
    theme(plot.title = element_text(
      hjust=0.5, size=14,face="bold"))
  ggsave(file = paste0(name, "_UMAP_by_seuratCluster.png"), 
         dpi = 300, width = 7, height = 7)
  
  UMAPPlot(object, label = F) + 
  labs(title="Unsupervised Clustering") + 
    theme(plot.title = element_text(
      hjust=0.5, size=14,face="bold"))
  
  ggsave(file = paste0(name, "_UMAP_by_seuratCluster_legend.png"), 
         dpi = 300, width = 7.7, height = 7)
  
  if (MyIdent1 != "") {
    Idents(object) <- MyIdent1
    num = Idents(object) %>% unique %>% length()
    
    UMAPPlot(object, label = F)
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, ".png"), 
           dpi = 300, width = 7, height = 7)
    
    UMAPPlot(object, group.by = MyIdent1, label = F)
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, ".png"), 
           dpi = 300, width = 7, height = 7)
    
    UMAPPlot(object, group.by = MyIdent1, split.by = MyIdent1, 
             label = F) + NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent1, "_split.png"), 
           dpi = 300, width = 5 * num, height = 7)
    
    Idents(object) <- "seurat_clusters"
    UMAPPlot(object, label = F, split.by = MyIdent1) + NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster_split_by_", 
                         MyIdent1, ".png"), dpi = 300, width = 5 * num, height = 7)
  }
  
  if (MyIdent2 != "") {
    Idents(object) <- MyIdent2
    num = Idents(object) %>% unique %>% length()
    
    UMAPPlot(object)
    ggsave(file = paste0(name, "_UMAP_", MyIdent2, ".png"), 
           dpi = 300, width = 7, height = 5)
    
    
    UMAPPlot(object, split.by = MyIdent2)
    ggsave(file = paste0(name, "_UMAP_by_split_by_", MyIdent2, 
                         ".png"), dpi = 300, width = 5 * num, height = 5)
    UMAPPlot(object, group.by = MyIdent2)
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent2, ".png"), 
           dpi = 300, width = 7.7, height = 7)
    UMAPPlot(object, group.by = MyIdent2, split.by = MyIdent2) + 
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_by_", MyIdent2, "_split.png"), 
           dpi = 300, width = 5 * num, height = 5)
  }
}




####################################
## Cluster Annotation後の check graph
## Usage: 
## fun.Annotated.Clustering.Check.Graph(ECMC, "Annotation1", "State", "ECMC_dim3Res0.5")
fun.Annotated.Clustering.Check.Graph <- 
function(object, MyIdent1, MyIdent2, samplename){
    Idents(object) <- MyIdent1
    num <- nlevels(object)
    name <- paste0("ClusterAnnotated_", samplename)
    
    #1) vlnplot
    VlnPlot(object, 
            c("nFeature_RNA", "nCount_RNA", "percent.mt"), #or percent.mt
            ncol=1, pt.size=0.1)
    ggsave(file = paste0(name, "_CountFeature_Vlnplot.png"),
           dpi = 300, width = 5, height = 15)
    
    # 3) UMAPで確認
    # ident: 自分でannotateしたやつ
    # ラベルなし，legendあり
    UMAPPlot(object, label=FALSE) +
    labs(title="Cluster Annotation") + 
    theme(plot.title = element_text(
          hjust=0.5, size=14,face="bold"))
    ggsave(file = paste0(name, "_UMAP_nolabel.png"),
      dpi = 300, width = 8, height = 7) 
    
    # ラベルあり，legened　なし
    UMAPPlot(object, label=TRUE) + 
      labs(title="Cluster Annotation") + 
      theme(plot.title =element_text(
            hjust=0.5, size=14, face="bold")) +
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_nolegend.png"),
           dpi = 300, width = 7, height = 7) 
    
    # annotateしたclusterを他の任意のclusterでsplitする
    # ラベルなし．legend あり
    if(MyIdent2 !=""){
    Idents(object) <- MyIdent2
    num2 = Idents(object) %>% unique() %>% length()
    
    Idents(object) <- MyIdent1
    UMAPPlot(
      object, 
      split.by=MyIdent2, 
      label=FALSE) +
      theme(
        plot.title = 
          element_text(hjust=0.5, size=14, face="bold"))
    ggsave(file = paste0(name, "_UMAP_by_seuratCluster_split_", MyIdent2, ".png"),
           dpi = 300, width = 5*num2, height = 5) 
    }
}




####################################
## Cluster Annotation後の check graph
## Usage: 
## fun.Annotated.Clustering.Check.Graph(ECMC, "Annotation1", "State", "ECMC_dim3Res0.5")
fun.Annotated.Clustering.Check.Graph.v2 <- 
  function(object, MyIdent1, MyIdent2, samplename){
    
    Idents(object) <- MyIdent1
    num <- nlevels(object)
    name <- paste0("ClusterAnnotated_", samplename)
    
    #1) vlnplot
    VlnPlot(object, c("nFeature_RNA", "nCount_RNA", 
                      "percent.mt", "percent.ribo",
                      "S.Score", "G2M.Score"), 
            ncol = 1, pt.size = 0.1) +
      labs(x="") 
    ggsave(file = paste0(name, "_CountFeature_Vlnplot.png"), 
           dpi = 300, width = 4, height = 16)
    
    #2) FeaturePlot
    FeaturePlot(object, c("nFeature_RNA", "nCount_RNA", 
                          "percent.mt", "percent.ribo",
                          "S.Score", "G2M.Score"), ncol = 2)
    ggsave(file = paste0(name, "_CountFeature_Featureplot.png"), 
           dpi = 300, width = 8, height = 8)
    
    # 3) UMAPで確認
    # ident: 自分でannotateしたやつ
    # ラベルなし，legendあり
    UMAPPlot(object, label=FALSE) +
      labs(title="Cluster Annotation") + 
      theme(plot.title = element_text(
        hjust=0.5, size=14,face="bold"))
    ggsave( file = paste0(name, "_UMAP_nolabel.png"),
            dpi = 300, width = 8, height = 7) 
    
    # ラベルあり，legened　なし
    UMAPPlot(
      object,
      label=TRUE) + 
      labs(title="Cluster Annotation") + 
      theme(
        plot.title =
          element_text(
            hjust=0.5, size=14, face="bold")) +
      NoLegend()
    ggsave(file = paste0(name, "_UMAP_nolegend.png"),
           dpi = 300, width = 7, height = 7) 
    
    # annotateしたclusterを他の任意のclusterでsplitする
    # ラベルなし．legend あり
    if(MyIdent2 !=""){
      Idents(object) <- MyIdent2
      num2 = Idents(object) %>% unique() %>% length()
      
      Idents(object) <- MyIdent1
      UMAPPlot(
        object, 
        split.by=MyIdent2, 
        label=FALSE) +
        theme(
          plot.title = 
            element_text(hjust=0.5, size=14, face="bold"))
      ggsave(file = paste0(name, "_UMAP_by_seuratCluster_split_", MyIdent2, ".png"),
             dpi = 300, width = 5*num2, height = 5) 
    }
  }

