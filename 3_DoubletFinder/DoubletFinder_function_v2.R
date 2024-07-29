## DoubletFinder test
library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(fs)
library(DoubletFinder)
packageVersion(pkg = "Seurat")


## updated on 05/14/24

## https://rpubs.com/kenneditodd/doublet_finder_example
## PCの数を自動的に判定

## QC (low qualityの細胞を除去)後にこの関数を通す


### Normalize data -> Scale dataで作った場合．#####################################################
Fun.doublet.finder.normal.v2 <- function(serutat_obj, dimention, doublet_percentage){
  
  require(DoubletFinder)
  require(fs)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.serutat_obj <- paramSweep(serutat_obj, PCs = 1:dimention) #sctはなくなった？
  sweep.serutat_obj <- summarizeSweep(sweep.serutat_obj, GT = FALSE)
  bcmvn_serutat_object <- find.pK(sweep.serutat_obj)
  
  # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn_serutat_object[which.max(bcmvn_serutat_object$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    print(paste0("Optimal pk is set to be ", optimal.pk))
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- serutat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(doublet_percentage*nrow(serutat_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  serutat_object_dbl <- doubletFinder(
    serutat_obj, 
    PCs = 1:dimention, 
    pN = 0.25, 
    pK = optimal.pk, 
    nExp = nExp_poi.adj, 
    reuse.pANN = FALSE)#, 
  #sct = FALSE)
  
  pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  DimPlot(serutat_object_dbl, group.by=pANNname)
  ggsave("DoubletFinder_object_1stround.png", width=8, height=8, dpi = 300)
  #   serutat_object_dbl@meta.data$DF.classifications_0.25_0.09_701 %>% table()
  # Doublet Singlet 
  # 701    8647 
  
  # Seurat metadata column name for previously-generated pANN results. Argument should be set to FALSE (default) for initial DoubletFinder runs. Enables fast adjusting of doublet predictions for different nExp.
  # pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  # serutat_object_dbl <- doubletFinder(
  #   serutat_object_dbl, 
  #   PCs = 1:dimention, 
  #   pN = 0.25, 
  #   pK = optimal.pk, 
  #   nExp = nExp_poi.adj, 
  #   reuse.pANN = pANNname)#, 
  # #sct = FALSE)
  # 
  # pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  # DimPlot(serutat_object_dbl, group.by=pANNname)
  # ggsave("DoubletFinder_object_2ndround.png", width=8, height=8, dpi = 300)
  # #serutat_object_dbl@meta.data$DF.classifications_0.25_0.09_648 %>% table()
  # 
  return(serutat_object_dbl)
}
