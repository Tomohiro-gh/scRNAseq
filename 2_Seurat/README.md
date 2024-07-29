# Seurat
### Seuratのエラーなどメモしておく

## Install
```r
remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
```

## version管理
#### check version of seurat
```r
## check suerat version
packageVersion("Seurat")
## checke object version
Version(pbmc_small)
```
  
-------------
## Dataの読み込み
#### 複数サンプルを一括で一気に読み込む
参考： https://zenn.dev/rchiji/articles/f92da3eef917ab

```r
dir_rawdata　<- "path/to/countmatrix"
  ## ここではディレクトリ内のfiltered_feature_bc_matrix．h５にマッチするファイルを読み込んでいく
(paths <- dir(path = dir_rawdata, pattern = "filtered_feature_bc_matrix.h5"))


## Object listの作成
obj_list <- lapply(paths, function(x){
  
  # 細胞バーコードに付与するサンプル名の用意。ファイルパスの中からサンプル名の箇所を取得している。
  GSMname <- stringr::str_split(string = x, pattern = "_", simplify = T)[,1]
  sample_name <- stringr::str_split(string = x, pattern = "_", simplify = T)[,2]
  fname <- stringr::str_split(string = x, pattern = "_filtered_feature_bc_matrix.h5", simplify = T)[,1]
  
  ## 1文字目がYで始まれば young, OならOldと名前をつける
  if(substr(sample_name, 1, 1) == "Y"){
    age <- "Young"
    }
  if(substr(sample_name, 1, 1) == "O"){
    age <- "Old"
  }
  
    # fname <- gsub(pattern = "_rna.h5", replacement = "", x = fname)
  x <- Read10X_h5(paste0(dir_rawdata, x))
  
  x <- CreateSeuratObject(counts = x, project = fname)
  
  x <- AddMetaData(object = x, metadata = fname, col.name = "SampleName")
  x <- AddMetaData(object = x, metadata = age, col.name = "Age")
  x <- AddMetaData(object = x, metadata = "Kidney", col.name = "Tissue")
 
  x <- RenameCells(x, fname)
  
  return(x)
})

```

-------------

## ロードしたオブジェクトのエラー
```r
## Error: Please run UpdateSeuratObject on your object
## Error during wrapup: no slot of name "images" for this object of class "Seurat"
## Error: no more error handlers available (recursive errors?); invoking 'abort' restart
```
Solution:
Converting this matrix to Seurat object #4515
https://github.com/satijalab/seurat/issues/4515
```r
  seurat_obj <- UpdateSeuratObject(seurat_obj)
```


## metadataの追加：既存のmetaを参考にして新たに追加するversion
```r
ECatlasQC@meta.data <- 
  ECatlasQC@meta.data %>% 
  mutate(rough_annotation =
           case_when(Cluster %in% A ~ 'Artery',
                     Cluster %in% V ~ 'Vein',
                     Cluster %in% C ~ 'Capillary',
                     Cluster %in% M ~ 'Mitotic',
                     Cluster %in% An ~ 'Angiogenic',
                     Cluster %in% L ~ 'Lymphatic'))


#また，一部のみを変更し，以下は同様としたい場合は， # TRUE ~ で指定する
SeuratObj@metadata <-
  SeuratObj@metadata %>% 
  mutate(
    Modified_MCannotation =
      case_when(rownames(AnnotationData) %in% Barcode_SMC ~ 'SMC',
                TRUE ~ Author_Annotation))

```
