## Clusteringをどれくらいの数値で行うかざっくりみるプロット
# 1 : Elbow
# 2 : Heatmap
# 3: JackStraw plot

## Requirement 下記の３行が実施されていること
  #DefaultAssay(SeuratObject) <- "integrated"
  #SeuratObject <- ScaleData(SeuratObject, verbose=FALSE) #1
  #SeuratObject <- RunPCA(SeuratObject, npcs=75, verbose=FALSE) #2

## Fun.Elbow.Heat.JackStraw(SeuratObject,filename)
## Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/Dim_Determination.R


Fun.Elbow.Heat.JackStraw <- function(SeuratObject, filename){
  
  DefaultAssay(SeuratObject) <- "integrated"
  
# 1) elbow plot :  横軸が主成分の次元を表し、降下曲線から水平になる点（Elbow）が最適な主成分の次元数
elb <- ElbowPlot(SeuratObject, ndims = 50) + 
    labs(title=paste0("elbow plot : ", filename)) + 
    theme(plot.title = element_text(hjust = 0.5))

# 2)  heatmap
# 下記でかく
# PC20までを一度に可視化 -> 二分されているように見えなければ、その主成分にはもう情報がないと判断

# 3): Jack-Straw plot (計算量が非常に多い): 30分くらいかかってしまう
# 黒点線で示されたランダム状態より上に振れているほど情報量が多い
#SeuratObject <- JackStraw(SeuratObject, num.replicate = 100)
#SeuratObject <- ScoreJackStraw(SeuratObject, dims = 1:20)
#options(repr.plot.width = 8, repr.plot.height = 6)
#jsp <- JackStrawPlot(SeuratObject, dims = 1:20) + 
#  labs(title=paste0("JackStraw Plot PC1-40 : ", filename)) + 
#      theme(plot.title = element_text(hjust = 0.6))
    
pdfname <- paste0("DimentionDetermination-", filename, ".pdf")

pdf(pdfname, width = 10, height = 8)

#1 Elbow plot     
plot(elb)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)

#2 Dim Heatmap: 15こずつ，書いてく 
# Dim 1-20
DimHeatmap(SeuratObject, dims = 1:20, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = "filename", cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
# Dim 21-40
DimHeatmap(SeuratObject, dims = 21:40, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = filename, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
# Dim 41-60
DimHeatmap(SeuratObject, dims = 41:60, cells = 500, balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = filename, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)

#3 Jack Straw Plot
#plot(jsp)

dev.off()  

}


Fun.Elbow.Heat.JackStraw(ed, "Integration_2xUW_3xWO")
