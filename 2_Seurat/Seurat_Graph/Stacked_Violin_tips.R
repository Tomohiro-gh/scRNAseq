## Seurat カラーを取得

library(Seurat)
library(scales)
library(ggsci)

## get length
Idents(seurat_obj) %>% unique() %>% length() #19

## color に関して
## https://paulvanderlaken.com/2020/03/20/how-to-standardize-group-colors-in-data-visualizations-in-r/

show_col(hue_pal()(19))
## or 
# color をvectorへ入れる
myColors = hue_pal()(19)

## またはggsciのカラーコードを使う場合，　例えばaaas colorにしてみると
    library(ggsci)
    show_col(pal_aaas("default")(10))
myColors_ggsci = pal_aaas(palette = c("default"), alpha = 1)(10)


## Identの変更
Idents(seurat_obj) <- "CellAnnotation_Fine"


## stacked violin plotで色をクラスターごとにしたい場合
## split.byに Identを入れておくといい
## stacked violinについての記載: https://zhuanlan.zhihu.com/p/652377364

VlnPlot(object = seurat_obj, 
        features = MarkerGenes, 
        split.by = "CellAnnotation_Fine", 
        cols = myColors,
        stack = TRUE,
        flip = TRUE,
        pt.size = 0.1) + 
    labs(x = "") +
    myTHEME
