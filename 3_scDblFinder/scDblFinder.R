## scdbl finder analysis
library(Seurat)
library(SeuratData)
library(scDblFinder)
library(tidyverse)

  pbmc
  # An object of class Seurat 
  # 12627 features across 500 samples within 1 assay 
  # Active assay: RNA (12627 features, 1940 variable features)
  # 3 layers present: counts, data, scale.data
  # 2 dimensional reductions calculated: pca, umap
  
## scDBlFinder
## https://bioconductor.org/packages/release/bioc/html/scDblFinder.html
  
#> 2 scDBl Finder
## https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html
  
# sce$scDblFinder.score : the final doublet score
# sce$scDblFinder.class : the classification (doublet or singlet)
  
# Run scDbltFinder 
sce <- as.SingleCellExperiment(pbmc)
sce <- sce |>
  scDblFinder(clusters = "seurat_clusters",
              dbr.per1k = 0.004) # 10Xの場合，設定しなくてOK, 1000 = 0.8% に設定されている



### cluster関係なくdoubletを検出したい場合．
# we run scDblFinder (providing the unusually high doublet rate)
sce <- as.SingleCellExperiment(pbmc)
sce <- scDblFinder(sce, dbr = 0.1) # clusterは指定しない

  sce$scDblFinder.score %>% head
  table(sce$scDblFinder.class)
  sce@colData@listData %>% as.data.frame() %>% head()
  
  
  
# Explore results and add to seurat object
meta_scdblfinder <- 
  sce@colData@listData %>%
  as.data.frame() %>% 
  dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
    head(meta_scdblfinder)

rownames(meta_scdblfinder) <- sce@colData@rownames
  head(meta_scdblfinder)
pbmc <-
  pbmc %>% 
  AddMetaData(
    metadata = meta_scdblfinder %>% dplyr::select('scDblFinder.class'))

    head(pbmc@meta.data)
    table(pbmc$scDblFinder.class)
    
  pbmc %>% 
    DimPlot(group.by = "scDblFinder.class")    
  
  