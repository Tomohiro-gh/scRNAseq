library(ggplot2)
library(openxlsx)
library(stringr)

## Dot plot 自作 version
#example data
GOres <- read.xlsx("/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp130_WH_scRNAseq_Analysis/SeuratAnalysis_v5/2-4_Integration_SCT_thenRBCremoval/220601_GOKEGG_ClusterProfiler/3_StromalCell/ClusterProfiler_Stromal cell_vs2dpi.xlsx")

GOres$Condition <- "vs2dpi"


# GeneRatioは xx/yyyyの形になっているので使えない．．．
## / の前と後で文字を抽出．１行ずつ行う
ratio <- GOres$GeneRatio
# new ratio　
new_ratio <- NULL

for(i in 1:length(ratio)){

  bunshi <- as.double(strsplit(ratio, "/")[[i]][[1]])
  bunbo <- as.double(strsplit(ratio, "/")[[i]][[2]])
  new_ratio[i] <- bunshi/bunbo
}

GOres$new_ratio <- new_ratio



# plot: dot plot
dp <- ggplot(data = GOres,
             aes(x = Condition,
                 y = Description, 
                 color = `p.adjust`,
                 size = new_ratio)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(plot.margin= unit(c(1, 1, -1, 1), "lines")) + 
  ggtitle("GO enrichment analysis")

plot(dp)