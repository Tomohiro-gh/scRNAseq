library(Seurat)
library(sctransform)
library()


## Integration function

## CCA integration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Fun.Integration.CCA <- function(merged_obj, split_condition){
  # split the dataset into a list of two seurat objects (stim and CTRL)
  seurat_list <- SplitObject(merged_obj, split.by = split_condition)
  
  # normalize and identify variable features for each dataset independently
  seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = seurat_list)
  
  ## Perform Integration
  Anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
  # this command creates an 'integrated' data assay
  seurat_obj_int <- IntegrateData(anchorset = Anchors)
  
  return(seurat_obj_int)
  
}

  ## Perform integrated analysis
  DefaultAssay(seurat_obj_int) <- "integrated"
  seurat_obj_int <- ScaleData(seurat_obj_int, verbose = FALSE)
  seurat_obj_int <- RunPCA(seurat_obj_int, npcs = 30, verbose = FALSE)
  seurat_obj_int <- RunUMAP(seurat_obj_int, reduction = "pca", dims = 1:30)
  seurat_obj_int <- FindNeighbors(seurat_obj_int, reduction = "pca", dims = 1:30)
  seurat_obj_int <- FindClusters(seurat_obj_int, resolution = 0.5)
  d

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  
## RPCA integration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  
  
  

## SCtransform integration >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  