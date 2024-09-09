library(Seurat)
library(ggplot2)



# ## Feature plot function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fun.fp.batch.v3.1 <- function (object, genelist, Redec, samplename, row.num) {
  require(patchwork)
  
  my_theme <-  theme(text = element_text(hjust=0.5, size=20, face="bold"),
                     axis.text = element_text(hjust=0.5, size=16, face="bold"),
                     axis.ticks.length = unit(2, "mm"),
                     axis.line = element_line(colour = "black", linewidth=1.0),
                     aspect.ratio = 1)
  
  ps <- list()
  for (gene in genelist) {
    p <- FeaturePlot(object, gene, reduction = Redec, cols = c("#c0c0c0", "#dc143c"))
    ps[[gene]] <- p
  }
 # genenames <- paste(genelist, collapse = ",")
    col.num <- ceiling(length(genelist)/row.num)
    wid <- 4.5 * col.num
    hei <- 4.5 * row.num
    wrap_plots(ps, nrow = row.num)
  ggsave(paste0("2_ClusteringCondition_FP_", samplename, ".png"), 
         dpi = 300, width = wid, height = hei)
}
# ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
