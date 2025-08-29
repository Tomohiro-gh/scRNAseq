#> 08/05/25
library(DoubletFinder)
library(scDblFinder)
library(Seurat)
library(SeuratData)
library(tidyverse)
#> Tutorial
#> https://github.com/chris-mcginnis-ucsf/DoubletFinder
#> 
#> Scdble finder で doublet判定 > Doublet finderの流れ



#> とにかく最初にcluster annotationまで行う

InstallData('pbmc3k') # run only at first time
pbmc <- LoadData('pbmc')
  pbmc
  # An object of class Seurat 
  # 12627 features across 500 samples within 1 assay 
  # Active assay: RNA (12627 features, 1940 variable features)
  # 3 layers present: counts, data, scale.data
  # 2 dimensional reductions calculated: pca, umap
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  VlnPlot(pbmc, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3)
  
pbmc <- subset(pbmc, 
               subset = nFeature_RNA > 200 &
                 nFeature_RNA < 2500 &
                 percent.mt < 5)

VlnPlot(pbmc, 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)


pbmc <- 
  pbmc |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |> 
  RunUMAP(dims = 1:10) |>
  FindNeighbors(dims = 1:10) |>
  FindClusters(resolution = 0.5)

  pbmc |>
    DimPlot(reduction = "umap", label = TRUE) # confirmation

# cluster annotation  
new.cluster.ids <-
    c("Naive CD4 T", 
      "CD14+ Mono", 
      "Memory CD4 T",
      "B",
      "CD8 T",
      "FCGR3A+ Mono",
      "NK", 
      "DC", 
      "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
pbmc$CellType <- Idents(pbmc)

  pbmc |>
    DimPlot(reduction = "umap", label = TRUE, pt.size = 0.5) +
    NoLegend()

  
#> SCdbl Finder ------------------
sce <- pbmc |>
    as.SingleCellExperiment() |>
    scDblFinder(#clusters = "seurat_clusters",
                samples = NULL,
                dbr.per1k = 0.004) 
  # 10Xの場合，設定しなくてOK, 1000 = 0.8% に設定されている
  
  #> confirmation
sce$scDblFinder.score %>% head
  table(sce$scDblFinder.class)
  sce@colData@listData %>% as.data.frame() %>% head()
  
#> join scdblfinder
meta_scdblfinder <- 
    sce@colData@listData %>%
    as.data.frame() %>% 
    dplyr::select(starts_with('scDblFinder')) # 'scDblFinder.class')
      head(meta_scdblfinder)
  
rownames(meta_scdblfinder) <- sce@colData@rownames
  head(meta_scdblfinder)

pbmc <-
  pbmc |> 
    AddMetaData(
      metadata = meta_scdblfinder %>% 
        dplyr::select('scDblFinder.class'))
  
  head(pbmc@meta.data)
  table(pbmc$scDblFinder.class)
  
#> ここでsingletとdoubletの表記を書き換えておく
pbmc@meta.data <- 
  pbmc@meta.data |>
  dplyr::mutate(
    scDblFinder.class = 
      str_c(str_to_upper(
        str_sub(scDblFinder.class, 1, 1)),
        str_sub(scDblFinder.class, 2, -1)))
    table(pbmc$scDblFinder.class) #確認

  # confirmation
  pbmc %>% 
    DimPlot(group.by = "scDblFinder.class")  



mesage("ここからDoubletFinderを使った設定を始めるよ")
require(DoubletFinder)  
require(fs)


## pK Identification (no ground-truth) ----------------------
sweep.res.list <- paramSweep(pbmc, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

data(pbmc_small)
seu <- pbmc_small
sweep.list <- paramSweep(seu)
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

## 現在はうまくグラフが出ないので．自作してみる.




## pK Identification (ground-truth) -----------------------------
sweep.res <- pbmc |>
  paramSweep(PCs = 1:10, sct = FALSE)
gt.calls <- pbmc@meta.data[rownames(sweep.res[[1]]), "scDblFinder.class"]
## GT is a vector containing "Singlet" and "Doublet" calls recorded using sample multiplexing classification and/or in silico geneotyping results 
sweep.stats <- 
  sweep.res |>
  summarizeSweep(GT = TRUE,
                 GT.calls = gt.calls)

bcmvn_kidney <- find.pK(sweep.stats)

## Homotypic Doublet Proportion Estimate --------------------
annotations <- pbmc@meta.data$CellType
homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.075*nrow(pbmc@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1 - homotypic.prop))

## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
pbmc <- 
  pbmc |>
  doubletFinder(
    PCs = 1:10,
    pN = 0.25,
    pK = 0.09,
    nExp = nExp_poi,
    # reuse.pANN = NULL,
    sct = FALSE)

Ln <- length(colnames(pbmc@meta.data))
pANN_name <- colnames(pbmc@meta.data)[Ln]

#> 250806 この段階でエラー　おそらくバグ
# pbmc <- 
#   pbmc |>
#   doubletFinder(
#     PCs = 1:10,
#     pN = 0.25,
#     pK = 0.09, 
#     nExp = nExp_poi.adj,
#     reuse.pANN = pANN_name,
#     annotations = annotations,
#     sct = FALSE)
  
#> manually worked 
#>https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/220#issuecomment-2703992695

pN <- 0.25
pK <- 0.01
nExp <- nExp_poi.adj
reuse.pANN <- pANN_name
pANN.old <- pbmc@meta.data[ , reuse.pANN]
classifications <- rep("Singlet", length(pANN.old))
classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
pbmc@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications

DimPlot(pbmc, group.by = "DF.classifications_0.25_0.09_198")
DimPlot(pbmc, group.by = "DF.classifications_0.25_0.01_164")　# 学習あり


bcmvn.max <- 
  sweep.serutat_obj[which.max(sweep.serutat_obj$BCmetric), ]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  message(paste0("Optimal pk is set to be ", optimal.pk))
  

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
