library(Seurat)
library(wesanderson)
library(patchwork)


message("loading function: 
  FUN.FeatureP.batch.MarkerGene(
      seurat_obj, genelist, Redec, samplename, row.num)")
### Feature plot function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
FUN.FeatureP.batch.MarkerGene <- function(seurat_obj, genelist, Redec, samplename, row.num) {
  require(patchwork)
  
  my_theme <-  theme(text = element_text(hjust=0.5, size=20, face="bold"),
                     axis.text = element_text(hjust=0.5, size=16, face="bold"),
                     axis.ticks.length = unit(2, "mm"),
                     axis.line = element_line(colour = "black", linewidth=1.0),
                     aspect.ratio = 1)
  ## gradient color
  my_pal = wes_palette("Zissou1", 20, type = "continuous")
  
  ps <- list()
  
  for (gene in genelist) {
    p <- FeaturePlot(seurat_obj, gene, reduction = Redec) + 
      scale_color_gradientn(colours = my_pal) +
      my_theme
    
    ps[[gene]] <- p
  }
  
  # genenames <- paste(genelist, collapse = ",")
  col.num <- ceiling(length(genelist)/row.num)
  wid <- 4.415 * col.num
  hei <- 4.415 * row.num
  
  wrap_plots(ps, nrow = row.num)
  
    ggsave(paste0("2_ClusteringCondition_", samplename, ".png"), 
            dpi = 300, width = wid, height = hei)
}
### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
###

#> FUN.FeatureP.batch.MarkerGene(seurat_obj, genelist, Redec, samplename, row.num)