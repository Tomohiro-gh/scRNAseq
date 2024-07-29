
## Setting Idents
# sanoke code: 210130_DiffNN_fromSourceFile.R
# Ident系の操作
celltype_id = "State_Annotation1" # metadata column name of the cell type of interest
seurat_obj = 
  SetIdent(
    seurat_obj,
    value = seurat_obj[[celltype_id]])

## 
seurat_obj[[]] %>% head()
## Stateに UW or WOが入っている

## ２つに分ける
Idents(seurat_obj) <- "State"

UW = seurat_obj %>% subset(idents = "UW")
WO = seurat_obj %>% subset(idents = "WO")

## それぞれで 解析




## Identの確認
seurat_obj[[]] %>% head()
seurat_obj[[]] %>% colnames
unique(Idents(seurat_obj))

## Levelの変更: 
seurat_obj$CellAnnotation_Fine <- 
  factor(x= seurat_obj$CellAnnotation_Fine,
         levels=c("Basal Epithelia", "Spinous Epithelia", "Mitotic Epithelia",
                  "Endothelial Cell", "Lymphatic EC", 
                  "Immune Cell",  "Macrophage", "Neutrophil", "T Cell",
                  "Fibroblast", "Pericyte", "SMC",  
                  "HF", "HFSC", "Subaceous Gland",
                  "Fast Skeltal Muscle","Muscle Satellite Cell",
                  "Shcwann Cell",
                  "EpiECMC?"))
seurat_obj@meta.data$CellAnnotation_Fine %>% unique() 