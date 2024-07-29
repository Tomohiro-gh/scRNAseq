
## Version1
DotPlot(obj, features = markers, cols=c("#5F4B8BFF", "#ED2B33FF"), assay = "RNA", col.min = 0.3, col.max = 0.8, dot.min=0.12, dot.scale = 1,
        cluster.idents=F)+
c+ 
  #scale_color_viridis_c(name = 'log2 (count + 1)') + 
  cowplot::theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 8, family="TT Times New Roman"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1, size = 8.5, family="TT Times New Roman"),
        legend.text = element_text(size=8),
        legend.title = element_text(size = 9)) + 
  scale_color_gradientn(colours = viridis::magma(20), limits = c(0,1), oob = scales::squish, name = 'log2 (count + 1)')


## Version2
## average Expressionの上限をmaxにする場合，
d <- DotPlot(obj, features = GOIs, group.by = "Annotation")


Obj %>% 
  DotPlot(., features =  GOIs, group.by = "Annotation") +
  scale_color_gradientn(colours = viridis::magma(20), 
                        limits = c(0, max(d$data$avg.exp)),
                        oob = scales::squish, 
                        name = 'log2 (count + 1)') +
  theme(axis.text.x = 
          element_text(angle = 45, vjust = 0.5, hjust=1, size = 20, family="TT Times New Roman"),
        axis.text.y = 
          element_text(angle = 0, vjust = 0.5, hjust=1, size = 20, family="TT Times New Roman"),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 9)) + 
  labs(x="", y="") + 
  RotatedAxis()
  #save
    ggsave("DotPlot.png", width = 8, height = 6, dpi=300)
