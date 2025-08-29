library(UCell)
library(Seurat)


HeartEC <- 
  HeartEC |> 
  UCell::AddModuleScore_UCell(
    features = list(RegAng = Ang),
    ncores = 14)
HeartEC |> 
  VlnPlot(features = "RegAng_UCell", 
          group.by = "ECannot",
          split.by = "Stage",
          cols  = color_stage,
          pt.size = 0) +
  geom_boxplot(position = position_dodge(width = 0.9), width=0.4) +
  # NoLegend() +
  labs(x="", y= "Score") +
  theme(text = element_text(size= 16, family = "Arial"),
        axis.text = element_text(size=20, family = "Arial", color = "black"),
        axis.text.y = element_text(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.ticks.length = unit(2, "mm"),
        axis.line = element_line(colour = "black", linewidth=0.5),
        aspect.ratio = 1.0)
ggsave("RegAnigogenesis_ECannot_UCell.pdf",
       width = 8, height = 7, dpi = 300)
