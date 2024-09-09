library(Seurat)
library(ggplot2)
library(cowplot)

## 発現を見るdotplot
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
FUN.DotPlot.Seurat.v1 <- function(object, Genelist, MyIdents, FileName = NULL){
  
  MyTheme <-  theme(
    axis.title.x = element_text(hjust = 0.5, size= 16, face = "bold", color = "black"),
    axis.title.y = element_text(hjust = 0.5, size= 16, face="bold", color = "black"),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 16, color = "black", family="TT Times New Roman"),
    axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, size = 16, color = "black", family="TT Times New Roman"),
    axis.line=element_line(colour = "black", linewidth=1.0),
    axis.ticks=element_line(colour = "black", linewidth = 1.0),
    axis.ticks.length = unit(2, "mm"),
    legend.text = element_text(size= 12),
    legend.title = element_text(size = 12)) 
  
  
  Idents(object) <- MyIdents
  
 DotPlot(object, 
          features = Genelist,
          assay = "RNA",
          scale = FALSE,
          #col.min = 0.3,
          #col.max = 0.8,
          #dot.min = 0.12,
          # dot.scale = 1,
          # cluster.idents=F
  ) + 
  # RotatedAxis() +
   coord_flip() +
    #scale_size(range = c(0, 5))+ 
    #scale_size_area(max_size = 5)+ 
    #scale_color_viridis_c(name = 'log2 (count + 1)') + 
    cowplot::theme_cowplot() +
    labs(x = "", y = "") +
    scale_color_gradientn(colours = viridis::magma(20)) + 
    MyTheme
    # scale_color_gradient2(low = "blue", high = "red", mid = "gray60")
    # scale_color_gradient2(low = "blue", high = "red", midpoint = mid)
  # scale_color_gradientn(colours = viridis::viridis(20))
  # scale_color_gradientn(colours = viridis::inferno(20)) #, #oob = scales::squish, name = 'log2 (count + 1)')
  
  
  ## Namingについて．Filenameの指定がなければ，Genelistの前から３つをfilenameにする
  if(is.null(FileName) == TRUE){
    if(length(Genelist) < 4){
      GeneNames <- paste(Genelist, collapse = ",")
    }else{
      GeneNames <- paste(Genelist[1:3], collapse = ",")
      FileName <- paste0(GeneNames, "_etc")
    }
  }
  
  ## 大きさの調整
  wid = 10 + length(unique(Idents(object))) * 0.2
  hei = 6 + length(Genelist) * 0.2
  ggsave(paste0("DotPlot_by_", MyIdents, "_", FileName, ".png"),
         width = wid, height = hei, dpi = 300, limitsize = FALSE)
}
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<





## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
##自作Dotplotのfunction : pvalueとFoldchangeで書いてみる
FUN.DotPlot.pval.vs.FC <- function(DEG_df, GeneList, colname_pval, colname_FC, colname_Gene, y_axis ,FileName = NULL){
  
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
  
  comp = DEG_df$Comparison %>% unique()
  #print(comp)
  df <- data.frame(Comparison = rep(comp, length(GeneList)),
                   Gene = GeneList)
  
  DEG_df <- DEG_df %>% 
    dplyr::rename(., p_val_adj = colname_pval) %>% 
    dplyr::rename(., log2FC = colname_FC) %>% 
    dplyr::rename(., Gene = colname_Gene) %>% 
    dplyr::rename(., Comparison = y_axis) %>% 
    dplyr::filter(Gene %in% GeneList) %>% 
    dplyr::left_join(df, DEG_df, by = c("Comparison", "Gene")) 
    dplyr::mutate(p_val_adj = 
                    ifelse(p_val_adj <= 1e-50, 1e-50,
                           p_val_adj))
  ## 並び順変更
  DFcombined$Comparison <- factor(x = DFcombined$Comparison,
                                  levels = comp)
  DFcombined$Gene <- factor(x = DFcombined$Gene,
                            levels = GeneList)  
  ## print(DFcombined$Comparison)
  ## 本体
  dp <- 
    ggplot(data = DFcombined,
           aes(x = Gene,
               y = Comparison)) + 
    # scale_size_continuous(limits = c(-3, 3)) +
    # scale_color_gradientn(colours = viridis::inferno(100),
    scale_color_gradientn(colours = viridis::viridis(100),
                          limits = c(min(DFcombined$log2FC), 
                                     max(DFcombined$log2FC))) + 
    scale_size_continuous(limits = c(0, 
                                     max(-log10(DFcombined$p_val_adj)))) +
    # plot.margin = unit(c(1, 1, -1, 1), "lines")) + 
    labs(title = "DEGs, pvalue x Fold Changes", x = "", y = "") + 
    theme_classic() +
    MyTheme
  
  dp <- ggplot(data = DFcombined,
               aes(x = Gene,
                   y = Comparison)) + 
    geom_point(data = subset(DFcombined, is.na(p_val_adj)),
               color = "white") + 
    geom_point(data = subset(DFcombined, !is.na(p_val_adj)),
               aes(color = avg_log2FC,
                   size = -log10(p_val_adj))) +
    
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
  
  ## Namingについて．Filenameの指定がなければ，Genelistの前から３つをfilenameにする
    if(is.null(FileName) == TRUE){
      if(length(Genelist) < 4){
        GeneNames <- paste(Genelist, collapse = ",")
      }else{
        GeneNames <- paste(Genelist[1:3], collapse = ",")
        FileName <- paste0(GeneNames, "_etc")
      }
    }
  
  ggsave(paste0("DotPlot_pval_x_FC_", FileName, ".png"), height = hei, width = wid, dpi = 300)
  
  }
