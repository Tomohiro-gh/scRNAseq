## Seurat でかくheatmap
library(Seurat)
library(ggplot2)
library()

## DoHeatmap
obj <- ScaleData(obj, features = rownames(obj))

DoHeatmap(object = obj,
          features = TopMarkers$gene,
          DoHeatmap(Object,
          group.by = MyIdent)  +
    scale_fill_gradientn(colors = viridis(100))
  )

    ggsave(paste0("DoHeatmap_Topmarkers_", n_genes, "_", FileName, ".png"),
              width = 20, height = 10, dpi = 300)
    
    
## dittoHeatmap
library(dittoSeq)

dittoHeatmap(obj,
             genes = TopMarkers$gene,
             order.by = "EC_int_cca_sub",
             scaled.to.max = FALSE,
             scale = "row")


    pheatmap::pheatmap(obj3[["RNA"]]$scale.data)
  
    
## 
library(ComplexHeatmap)
df <- obj[["RNA"]]$scale.data %>%
  as.data.frame() %>% 

mat <- obj[["RNA"]]$scale.data %>%
  as.matrix() %>% 
  colSort(obj$EC_int_cca_sub)


df <- df[DEGdf3_top20$Gene,]
  dim(df)
  
  order(colnames(obj$EC_int_cca_sub), levels(obj$EC_int_cca_sub))
  sort((obj$EC_int_cca_sub))
  
Heatmap(df, 
        name = "scale.data", #title of legend
        column_title = "Cells",
        row_title = "Features",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_order = obj$EC_int_cca_sub,
        row_order = DEGdf3_top20$Gene,,
        row_names_gp = gpar(fontsize = 7) # Text size for row names
)