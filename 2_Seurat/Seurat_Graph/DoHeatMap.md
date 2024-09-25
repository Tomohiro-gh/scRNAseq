
## DoHeatmap

##### update on 09/25/24
```r
## assume that draw heatmap of aggregate/average seurat object
FUN.DoHeatMap.v091324 <- function(SeuratAggr, DEG_features, slot){
  
  require(RColorBrewer)
  graphtitle = paste0( i, " reuglated genes in cluster ", j, ": ", k," slot")
 
  theme_heat <- theme(
    plot.margin= unit(c(0, 0, 0, 0), "inch"),
    text = element_text(hjust = 0.5, size = 16, face="bold", color = "black"),
    # axis.text.x = element_text(hjust=0.5, size=18, face="bold", color = "black"),
    axis.text.y = element_text(hjust=0.5, size=10, face="bold", color = "black")) 

  if(slot == "data"){
    color = viridis::magma(20)
  }else if(slot == "scale.data"){
    color = rev(brewer.pal(9, "RdBu"))
  }

  ## Heatmap
  DoHeatmap(SeuratAggr, 
            features = DEG_features,
            draw.lines = FALSE,
            slot = slot) +
    scale_fill_gradientn(colours = color) +
    labs(title = graphtitle) + 
    guides(colour=FALSE) + 
    theme_heat
  
  ggsave(paste0("DoHeart_", j, "_", i, "_regulated_",k ,"slot.png"),
         width = 10, height = 12, dpi = 300)
  }

```


--- 
#### How can I remove the legend for cell identity but not expression color bar using Doheatmap? #4557
 https://github.com/satijalab/seurat/issues/4557

 
- `DoHeatmap(...) + guides(colour=FALSE)`　で identを消す
- `DoHeatmap(...) + guides(fill=FALSE)` で　expressionを消す
- `DoHeatmap(...) + theme(legend.position="none")`　で両方のlegendをなくす

DoHeatmap(ss, features = m1$gene)+
  guides(colour=FALSE)

DoHeatmap(ss, features = m1$gene)+
  guides(fill=FALSE)
