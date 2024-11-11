library(Seurat)
library(ggplot2)
library(scales)

message("
 VlnPlot.Stacked.v1(
    seurat_object, genelist, MyIdent, FileName = NULL) ")
## 2 
VlnPlot.Stacked.v1 <- function(seurat_object, genelist, MyIdent, FileName = NULL){
  MyTheme <-  theme(# text = element_text(hjust=0.5, size=16, face="bold"),
    axis.text = element_text(hjust=0.5, size= 16, face = "bold"),
    axis.text.x = element_text(size = 20, angle = 60, vjust = 1.0, color = "black"),
    axis.text.y = element_text(size= 16 , face="bold"),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    axis.ticks.length = unit(2, "mm"),
    axis.line = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none",
    legend.justification = "left",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(-10,5,0,0),
    panel.spacing = unit(0, "lines"),
    panel.background = element_blank(),
    panel.border = element_blank(),
    plot.background = element_blank(),
    plot.margin = margin(14, 7, 7, 7, "pt"),
  )
  
  
  Idents(seurat_object) <- MyIdent
  myColors <- scales::hue_pal()(length(unique(Idents(seurat_object))))
  
  VlnPlot(object = seurat_object, 
          features = genelist,
          #split.by = MyIdent, 
          cols = myColors, 
          fill.by = "ident",
          stack = TRUE,
          flip = TRUE) + 
    labs(x = "") +
    MyTheme + 
    NoLegend()
  
  hei = 8 + length(genelist) * 0.4
  wid = 10 + length(unique(Idents(seurat_object))) * 0.2
  
  ## 名前
  if(length(genelist) >= 3){
    Genenames <- paste(genelist[1:3], collapse = ",")
    Genenames <- paste0(Genenames, ",etc")
  }else{Genenames = paste(genelist, collapse = ",")}
  
  imagename <- paste0("VlnStacked_", FileName,"_", Genenames, ".png")
  ggsave(imagename, dpi = 300, width = wid, height = hei,
         limitsize = FALSE)
  
}




# ## 
# FUN.Vln.Stacked.cca.cluster <- 
#   function (object, genelist, Filename) {
#     
#     MyTheme <-  theme(text = element_text(hjust=0.5, size=16, face="bold"),
#                       axis.text = element_text(hjust=0.5, size= 16, face="bold"),
#                       axis.text.x = element_text(angle = 90, vjust = 0.5, color = "black"),
#                       axis.title.y = element_blank(),
#                       axis.ticks.y = element_blank(),
#                       axis.text.y = element_blank(),
#                       axis.ticks.length = unit(2, "mm"),
#                       axis.line = element_line(colour = "black", linewidth=0.5),
#                       legend.position = "right",
#                       legend.justification = "left",
#                       legend.margin = margin(0,0,0,0),
#                       legend.box.margin = margin(-10,5,0,0),
#                       panel.spacing = unit(0, "lines"),
#                       panel.background = element_blank(),
#                       panel.border = element_blank(),
#                       plot.background = element_blank(),
#                       plot.margin = margin(14, 7, 7, 7, "pt"))
# 
#     kari_obj <- object
#     # myColors = scales::hue_pal()(length(unique(kari_obj$int_cca)))
#     kari_obj$int_cca <- factor(x=kari_obj$int_cca,
#                                levels=c("24", "12", 
#                                         "23", "6","10","15",
#                                         "2", "25","7","17", "3","14",
#                                         "5", "1", "0","4", "8", "19","20",
#                                         "11", "9","16", "18", "13", "21", "22"
#                                ))
#     
#     VlnPlot(object = kari_obj, 
#             features = genelist,
#             group.by = "int_cca",
#             # split.by = "CellType_AllCell_cca", 
#             # cols = myColors,
#             stack = TRUE,
#             flip = TRUE) + 
#       labs(x = "") +
#       MyTheme + 
#       NoLegend()
#     
#     wid = 6 + length(genelist)
#     ggsave(paste0("3_VlnStacked_", FileNameLabel,"_", Filename, ".png"), 
#            dpi = 300, width = wid, height = wid)
#     
#   }


