#> DotPlot
fun.DotPlot <- function(seurat_obj, genelist, Name = NULL, MyGroup){
  
  theme_dp <- theme(axis.text.y = element_text(size = 14, lineheight = 0.6),
                    axis.text.x = element_text(angle = 45, hjust = 1, size = 14))
  DotPlot(seurat_obj, 
          genelist,
          group.by = MyGroup) + 
    #split.by = "YoungorAged",
    #scale = FALSE
    coord_flip() +
    labs(x= "", y="") +
    theme_dp
  
  if(length(genelist)<=20){
    hei = 10
  }else{
    hei = length(genelist) * 0.5
  }
  ggsave(paste0("DotPlot_", Name, ".png"),
         width = 20, height = hei, dpi = 300)
} # end function ------



HeartCap |>
  DotPlot(features = Caps_representitive,
          group.by = "ECannot_Fine_DotPlot", 
          scale = TRUE) +
  scale_color_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrRd")) +
  labs(x="", y="") +
  theme(
    axis.title = element_text(family = "Arial"),
    axis.text = element_text(family = "Arial", size = 16),
    legend.title = element_text(family = "Arial", size = 8),
    legend.text = element_text(family = "Arial", size = 8),
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(size = 16,
                               angle = 90, hjust = 1, vjust= 0.5,
                               color = "black"))




#> Average Expression した後のDoeHeatmap
obj |>
  AverageExpression(
    group.by = "ECannot_Fine_x_Age",
    return.seurat = TRUE) |>
  DoHeatmap(features = GOIs,
            draw.lines = FALSE) + 
  NoLegend() +
  scale_fill_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  coord_fixed() +
  theme(
    axis.text = element_text(size = 12))


