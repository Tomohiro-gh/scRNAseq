library(Seurat)
library(ggplot2)
library(dplyr)

## ２つの遺伝子を指定し，その発現の関係を scatter plotで示す関数．

## Scatter plotをかくfunction
Fun.Scatter.Plot.2.genes_v1 <- function(TwoGenes, smooth){
  mtx = seurat_obj@assays[["RNA"]]@data[TwoGenes, ] %>% 
    as.matrix() %>% 
    t() 
  df = as.data.frame(mtx)
  p <- ggplot(data = df, aes_string(x=TwoGenes[1], y=TwoGenes[2])) + geom_point()
  
  if(smooth == TRUE){
    p <- p + geom_smooth()
  }
  
  p <- p + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_test() + 
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title.x = element_text(hjust=0.5, size= 16, face="bold"),
      axis.title.y = element_text(hjust=0.5, size= 16, face="bold"),
      axis.text.x = element_text(hjust=0.5, colour = "black", size= 14, face="bold"),
      axis.text.y = element_text(hjust=0.5, colour = "black", size= 14, face="bold"),
      #axis.line = element_line(colour = "black", linewidth = 0.75),
      axis.ticks = element_line(colour = "black"),
      text = element_text(size = 12))
  #legend.position="none")
  
  plot(p)
  
  ggsave(paste0("ScatterPlot_DataSlot_",TwoGenes[1], "_", y=TwoGenes[2],".png"),
         width=6, height=6, dpi=300)
}





## Scatter plotをかくfunction
Fun.Scatter.Plot.2.genes_v2 <- function(TwoGenes, smooth){
  mtx = seurat_obj@assays[["RNA"]]@scale.data[TwoGenes, ] %>% 
    as.matrix() %>% 
    t() 
  df = as.data.frame(mtx)
  p <- ggplot(data = df, aes_string(x=TwoGenes[1], y=TwoGenes[2])) + geom_point()
  
  if(smooth == TRUE){
    p <- p + geom_smooth()
  }
  
  p <- p + 
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_test() + 
    theme(
      plot.title = element_text(hjust=0.5, size=14, face="bold"),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title.x = element_text(hjust=0.5, size= 16, face="bold"),
      axis.title.y = element_text(hjust=0.5, size= 16, face="bold"),
      axis.text.x = element_text(hjust=0.5, colour = "black", size= 14, face="bold"),
      axis.text.y = element_text(hjust=0.5, colour = "black", size= 14, face="bold"),
      #axis.line = element_line(colour = "black", linewidth = 0.75),
      axis.ticks = element_line(colour = "black"),
      text = element_text(size = 12))
  #legend.position="none")
  
  plot(p)
  
  ggsave(paste0("ScatterPlot_ScaleDataSlot_",TwoGenes[1], "_", y=TwoGenes[2],".png"),
         width=6, height=6, dpi=300)
}
Fun.Scatter.Plot.2.genes_v2(TwoGenes = c("Cspg4", "Lgals1"), smooth = TRUE)
Fun.Scatter.Plot.2.genes_v2(TwoGenes = c("Cspg4", "Islr"), smooth = TRUE)

