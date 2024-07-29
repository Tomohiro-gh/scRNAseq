## DoubletFinder test

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(fs)
library(DoubletFinder)
packageVersion(pkg = "Seurat")



## https://github.com/chris-mcginnis-ucsf/DoubletFinder

# test_obj <- obj_list[[1]]
# 
# test_obj <- NormalizeData(test_obj)
# test_obj <- FindVariableFeatures(test_obj)
# test_obj <- ScaleData(test_obj)
# 
# test_obj <- RunPCA(test_obj)
# # ElbowPlot(obj, ndims = 50)
# # test_obj <- FindNeighbors(test_obj, dims = 1:30, reduction = "pca")
# # test_obj <- FindClusters(test_obj, resolution = 0.5, cluster.name = "seurat_clusters")
# test_obj <- RunUMAP(test_obj, dims = 1:30)
# 
# test_obj_DF <- Fun.doublet.finder.normal(test_obj, dimention = 30, doublet_percentage = 0.075)
# 
# dimention = 30
# sweep.test_obj <- paramSweep(test_obj, PCs = 1:30) #sctはなくなった？


### Normalize data -> Scale dataで作った場合．#####################################################
Fun.doublet.finder.normal <- function(serutat_obj, dimention, doublet_percentage){
  
  require(DoubletFinder)
  require(fs)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.serutat_obj <- paramSweep(serutat_obj, PCs = 1:dimention) #sctはなくなった？
  sweep.serutat_obj <- summarizeSweep(sweep.serutat_obj, GT = FALSE)
  bcmvn_serutat_object <- find.pK(sweep.serutat_obj)
  
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
    pK = 0.09, 
    nExp = nExp_poi, 
    reuse.pANN = FALSE)#, 
    #sct = FALSE)
  
    pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
      DimPlot(serutat_object_dbl, group.by=pANNname)
        ggsave("DoubletFinder_object_1stround.png", width=8, height=8, dpi = 300)
    #   serutat_object_dbl@meta.data$DF.classifications_0.25_0.09_701 %>% table()
      # Doublet Singlet 
      # 701    8647 
  
  # Seurat metadata column name for previously-generated pANN results. Argument should be set to FALSE (default) for initial DoubletFinder runs. Enables fast adjusting of doublet predictions for different nExp.
    pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  serutat_object_dbl <- doubletFinder(
    serutat_object_dbl, 
    PCs = 1:dimention, 
    pN = 0.25, 
    pK = 0.09, 
    nExp = nExp_poi.adj, 
    reuse.pANN = pANNname)#, 
    #sct = FALSE)
    
    pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
      DimPlot(serutat_object_dbl, group.by=pANNname)
        ggsave("DoubletFinder_object_2ndround.png", width=8, height=8, dpi = 300)
  #serutat_object_dbl@meta.data$DF.classifications_0.25_0.09_648 %>% table()
  
  return(serutat_object_dbl)
}





########################################################################
## ここまでは条件を決めてdimを決定しておく
serutat_obj <- SCTransform(serutat_obj)
serutat_obj <- RunPCA(serutat_obj)
ElbowPlot(serutat_obj, ndims = 50)
serutat_obj <- FindNeighbors(serutat_obj, dims = 1:40, reduction = "pca")
serutat_obj <- FindClusters(serutat_obj, resolution = 1, cluster.name = "sct_clusters")
serutat_obj <- RunUMAP(serutat_obj, dims = 1:40)
####. SCT version
Fun.doublet.finder.SCT <- function(serutat_obj, dims, doublet_percentage){
  
  require(DoubletFinder)
  require(fs)
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.serutat_obj <- paramSweep(serutat_obj, PCs = 1:dims, sct = TRUE)
  sweep.serutat_obj <- summarizeSweep(sweep.serutat_obj, GT = FALSE)
  bcmvn_serutat_object <- find.pK(sweep.serutat_obj)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  annotations <- serutat_obj@meta.data$sct_clusters
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  
  nExp_poi <- round(doublet_percentage*nrow(serutat_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  serutat_object_dbl <- doubletFinder(
    serutat_obj, 
    PCs = 1:40, 
    pN = 0.25, 
    pK = 0.09, 
    nExp = nExp_poi, 
    reuse.pANN = FALSE, 
    sct = TRUE)
  
  pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  DimPlot(serutat_object_dbl, group.by=pANNname)
  ggsave("DoubletFinder_object_1stround.png", width=8, height=8, dpi = 300)
  #   serutat_object_dbl@meta.data$DF.classifications_0.25_0.09_701 %>% table()
  # Doublet Singlet 
  # 701    8647 
  
  # Seurat metadata column name for previously-generated pANN results. Argument should be set to FALSE (default) for initial DoubletFinder runs. Enables fast adjusting of doublet predictions for different nExp.
  pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  serutat_object_dbl <- doubletFinder(
    serutat_object_dbl, 
    PCs = 1:40, 
    pN = 0.25, 
    pK = 0.09, 
    nExp = nExp_poi.adj, 
    reuse.pANN = pANNname, 
    sct = TRUE)
  
  pANNname = serutat_object_dbl@meta.data %>% colnames %>% tail(n=1)
  DimPlot(serutat_object_dbl, group.by=pANNname)
  ggsave("DoubletFinder_object_2ndround.png", width=8, height=8, dpi = 300)
  #serutat_object_dbl@meta.data$DF.classifications_0.25_0.09_648 %>% table()
  
  return(serutat_object_dbl)
}

