library(Seurat)


## Identの確認
seurat_obj[[]] %>% head()
seurat_obj[[]] %>% colnames
unique(Idents(seurat_obj))


# example code: 221223_CluaterAnnotation.R
#まずunsupervised clusterであるか確認
Idents(seurat_obj) %>% table()

#### Cluster annotation
## NEW cluster name
new.cluster.ids <- c("Fibroblast",               #0
                     "Basal",               #1
                     "Spinous",        #2
                     "Fibroblast",               #3
                     "ImmuneCell",               #4
                     "ImmuneCell",               #5
                     "Basal",               #6
                     "Basal",      #7
                     "Fibroblast",               #8
                     "ImmuneCell",   #9
                     "Mural",        #10
                     "HF",    #11
                     "Basal",               #12
                     "Prolif.Basal",    #13
                     "SubaceousGland",               #14
                     "Endothelial",        #15
                     "Unknown",        #16
                     "ImmuneCell",       #17
                     "Lympth ?",      #18
                     "SchwannCell",      #19
                     "ImmuneCell",            #20
                     "TCell") #21 

names(new.cluster.ids) <- levels(seurat_obj) 
seurat_obj <- 
  RenameIdents(
    seurat_obj, 
    new.cluster.ids)

#名前を変更した metadataをnew.meta.data.nameのmeta.dataへ格納する
seurat_obj$new.meta.data.name <- Idents(seurat_obj)



new.cluster.ids <- c("",   #0
                     "",               #1
                     "",        #2
                     "",               #3
                     "",               #4
                     "",               #5
                     "",               #6
                     "",      #7
                     "",               #8
                     "",   #9
                     "",        #10
                     "",    #11
                     "",               #12
                     "",    #13
                     "",               #14
                     "",        #15
                     "",        #16
                     "",       #17
                     "",      #18
                     "",      #19
                     "",            #20
                     "",        #21 
                     ""     #22
)


## Template for Cluster annotation
Idents(seurat_obj) <- "seurat_clusters" # unsupervised clustering後のクラスタレベル
nlevels(Idents(seurat_obj))


new.cluster.ids <- 
  c("",   #0
    "",   #1
    "",   #2
    "",   #3
    "",   #4
    "",   #5
    "",   #6
    "",   #7
    "",   #8
    "",   #9
    "",   #10
    "",   #11
    "",   #12
    "",   #13
    "",   #14
    "",   #15
    "",   #16
    "",   #17
    "",   #18
    "",   #19
    "",   #20
    "",   #21
    "",   #22
    "",   #23
    "",   #24
    "",   #25
    "",   #26
    "",   #27
    "",   #28
    "",   #29
    "",   #30
  ) 


names(new.cluster.ids) <- levels(seurat_obj) 
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# Author_Annotationというmeta.dataへidentを保存
seurat_obj$CellAnnotation <- Idents(seurat_obj)



# Identを変更
celltype_id = "CellAnnotation" # metadata column name of the cell type of interest
seurat_obj = 
  SetIdent(
    seurat_obj,
    value = seurat_obj[[celltype_id]])


## おまけ
table(
  seurat_obj@meta.data$Author_Annotation,
  seurat_obj@meta.data$age)



## levelの修正
Idents(seurat_obj) <-  "CellAnnotation"
levels(Idents(seurat_obj))

Idents(seurat_obj) <- 
  factor(Idents(seurat_obj),
         levels=c("", "", "",
                  "", "", "",
                  "", "", "",
                  "", "", "",
                  "", "", ""))

levels(Idents(seurat_obj))
