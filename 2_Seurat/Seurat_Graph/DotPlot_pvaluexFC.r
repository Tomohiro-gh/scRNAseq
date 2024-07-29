## 自作のdot plot: Fold changeとpvalueで表示

## data frameにfold changeと pvalue， sample名が入っていることが大事
FUN.DotPlot.Jisaku <- function(DEG_df, GeneList, FileName = NULL){
  
  DEG_df <- DEG_df %>% 
    filter(Gene %in% GeneList) %>% 
    dplyr::mutate(p_val_adj = 
                    ifelse(p_val_adj <= 1e-50, 1e-50,
                           p_val_adj))

  MyTheme <-  theme(
    axis.title.x = element_text(hjust = 0.5, size= 16, face = "bold", color = "black"),
    axis.title.y = element_text(hjust = 0.5, size= 16, face="bold", color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16, color = "black"),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 16, color = "black"),
    axis.line=element_line(colour = "black", linewidth=1.0),
    axis.ticks=element_line(colour = "black", linewidth = 1.0),
    axis.ticks.length = unit(2, "mm"),
    legend.text = element_text(size= 12),
    legend.title = element_text(size = 12)) #,
    
  dp <- ggplot(data = DEG_df,
               aes(x = Gene,
                   y = Comparison, 
                   color = avg_log2FC,
                   size = -log10(p_val_adj))) + 
    geom_point() + 
   # scale_size_continuous(limits = c(-3, 3)) +
   scale_color_gradientn(colours = viridis::inferno(100),
                          limits = c(min(DEG_df$avg_log2FC), 
                                     max(DEG_df$avg_log2FC))) + 
  scale_size_continuous(limits = c(0, 
                                   max(-log10(DEG_df$p_val_adj)))) +
      # plot.margin = unit(c(1, 1, -1, 1), "lines")) + 
    labs(title = "DEGs, pvalue x Fold Changes", x = "", y = "") + 
    theme_classic() +
    MyTheme
        
        plot(dp)
    
      hei_len = DEG_df$Comparison %>% unique() %>% length()
      wid = 10 + length(GeneList)*0.2
      hei = 6 + hei_len*0.2
      if(is.null(FileName) == TRUE){
        FileName = paste(GeneList, collapse = ",")
      }
        ggsave(paste0("DotPlot_pval_x_FC_", FileName, ".png"), height = hei, width = wid, dpi = 300)
}


## Eexcution
FUN.DotPlot.Jisaku(DEG_df = DEGL_list_ALL_combined_df, 
                   GeneList = ALLDEGs_UP_ligand, 
                   FileName = "UpregulatedGenes")
