library(Seurat)
library(SeuratWrappers)
library(slingshot)
library(cowplot)
library(ggplot2)
library(devtools)
library(dplyr)


## Dot plotで何色もgradientに色を使う

DotPlot(pbmc,
        features = markers,
        cols=c("#5F4B8BFF", "#ED2B33FF"),
        assay = "RNA", col.min = 0.3,
        col.max = 0.8,
        dot.min=0.12,
        dot.scale = 1,
        cluster.idents=F)+
  scale_size(range = c(0, 5))+ 
  scale_size_area(max_size = 5)+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 8, family="Arial"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 8.5, family="Arial"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) + 
  scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)')


################ Dot plot  ################################
fun.dp <- function(sub, markers, filename){
  pdfname <- paste0("DotplotEachCluster", filename, ".pdf")
  d <- DotPlot(sub, features=markers, cols=c("blue", "red","magenta", "cyan"), dot.scale=8, split.by="Condition") +
    RotatedAxis()
  pdf(pdfname, width=16, height=12)
  plot(d)
  dev.off()
}
# dotplot(wh.int.38.15, markers.dp, "Integration")


