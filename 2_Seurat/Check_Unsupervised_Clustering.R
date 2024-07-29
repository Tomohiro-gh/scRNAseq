# clustering後のcheck graph
# sampleのbatch effectの確認などに使用
### Graph: QC  
## Usage: fun.Clustering.Check.Graph(object, MyIdent1, MyIdent2, samplename)
## /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/Clustering_check_graph.R

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


# Example : 23/1/20
fun.Clustering.Check.Graph(ECPC6, "State", "cluster", "ECPC_dim6res0.5")

