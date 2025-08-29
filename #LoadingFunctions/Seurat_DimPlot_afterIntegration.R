print("fun.DimPlot.batch.v1 は resolutionを振った時に一度にplotするのに使用")

# ## Feature plot function >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
fun.DimPlot.batch.v1 <- function (object, 
                                  SeqList, 
                                  ReductionName = NULL,
                                  labelname=NULL, 
                                  row.num, 
                                  Filename) {
  require(patchwork)
  #> theme settings
  my_theme <- theme(text = element_text(hjust=0.5, size=16, face="bold"),
                    axis.text = element_text(hjust=0.5, size=16, face="bold"),
                    axis.ticks.length = unit(2, "mm"),
                    axis.line = element_line(colour = "black", linewidth=1.0),
                    aspect.ratio = 1)
  ##
  ps <- list()
  for (sl in SeqList) {
    labs_name = paste0(labelname, " ",  sl)
    p <- DimPlot(object, reduction = ReductionName, label = TRUE, group.by = sl) +
      labs(title = labs_name) + 
      my_theme
    ps[[sl]] <- p
  }
  
  col.num <- ceiling(length(SeqList)/row.num)
  wid <- 8 * col.num
  hei <- 8 * row.num
  wrap_plots(ps, nrow = row.num)
  ggsave(paste0("2_Dimplot_", Filename, ".png"), 
         dpi = 300, width = wid, height = hei)
}
# ## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
