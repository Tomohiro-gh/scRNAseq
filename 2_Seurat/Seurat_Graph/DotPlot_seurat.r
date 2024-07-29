## Sueuratのdot plot機能を使用．　colorをgradientにする

library(seurat)
library(ggplot2)
library(scales)

FUN.DotPlot.Seurat <- function(object, Genelist, MyIdents, FileName = NULL){
  
  # MyColor = c("#5F4B8BFF", "#ED2B33FF")
  # MyColor2 = c("#ffff00", "#4b0082")
  Idents(object) <- MyIdents
  
  DotPlot(object, 
          features = Genelist,
          #cols = MyColor,
          # group.by = MyIdents,
          assay = "RNA",
          scale = FALSE,
          #col.min = 0.3,
          #col.max = 0.8,
          #dot.min = 0.12,
          # dot.scale = 1,
          # cluster.idents=F
  ) + 
    RotatedAxis() + 
    #scale_size(range = c(0, 5))+ 
    #scale_size_area(max_size = 5)+ 
    #scale_color_viridis_c(name = 'log2 (count + 1)') + 
    cowplot::theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16, family="TT Times New Roman"),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 16, family="TT Times New Roman"),
          legend.text = element_text(size= 10),
          legend.title = element_text(size = 10)) + 
    labs(x = "", y = "") +
  scale_color_gradientn(colours = viridis::magma(20), oob = scales::squish, name = 'log2 (count + 1)')
  # scale_color_gradientn(colours = viridis::viridis(20), oob = scales::squish, name = 'log2 (count + 1)')
  # scale_color_gradientn(colours = viridis::inferno(20), oob = scales::squish, name = 'log2 (count + 1)')
  
  
  if(length(Genelist) < 4){
    GeneNames <- paste(Genelist, collapse = ",")
  }else{
    GeneNames <- paste(Genelist[1:3], collapse = ",")
    GeneNames <- paste0(GeneNames, "_etc")
  }
  wid = 10 + length(Genelist) * 0.2
  hei = 6 + length(unique(Idents(object))) * 0.2
  ggsave(paste0("DotPlot_by_", MyIdents, "_", FileName,"_", GeneNames, ".png"),
         width = wid, height = hei)
}
