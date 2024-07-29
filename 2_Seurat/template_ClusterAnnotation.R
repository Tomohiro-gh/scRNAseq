library(Seurat)


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
