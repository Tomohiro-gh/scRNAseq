

####### elbow plot
elb <- ElbowPlot(ECPC, ndims = 100) + 
  labs(title="ECPC standard integration") + 
  theme(plot.title = element_text(hjust = 0.5))
plot(elb)


### Dim heatmap
fun.Dim.Heat <- function(SeuratObject, dimention, Title){
  pdf("subclusterEC_DimHeat.pdf", width = 10, height = 8)
  DimHeatmap(
    SeuratObject, 
    dims = dimention, 
    cells = 500, 
    balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = Title, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
    dev.off()
}

fun.Dim.Heat(ECPC, 1:20, "ECPC, standard integration")


pdf("subclusterMC_DimHeat.pdf", width = 10, height = 8)
DimHeatmap(MC, dims = 1:20, cells = 500, balanced = TRUE)
mtext(side = 3, line=0.2, outer=T, text = "subclusterECMC", cex=1)
par(oma = c(0, 0, 1.5, 0)) 
options(repr.plot.width = 15, repr.plot.height = 15)
dev.off() 

pdf("subclusterECMC_DimHeat.pdf", width = 10, height = 8)
DimHeatmap(ECMC, dims = 1:20, cells = 500, balanced = TRUE)
mtext(side = 3, line=0.2, outer=T, text = "subclusterECMC", cex=1)
par(oma = c(0, 0, 1.5, 0)) 
options(repr.plot.width = 15, repr.plot.height = 15)
dev.off()  




## Dim, Resを変えながら，一気にclustering
fun.cluster.finding.UMAP <- 
  function(sub, filename){
  pdfname <- 
    paste0("ClusteringConditionUMAP_", filename, ".pdf")
  
  pdf(pdfname, width = 6, height = 6)
  
  for (i in 1:length(con_dim)){
    sub <- RunUMAP(sub, dims = 1:con_dim[i]) #3
    sub <- FindNeighbors(object = sub, dims = 1:con_dim[i]) #4
    for (j in 1:length(con_res)){
      sub <- FindClusters(object = sub, resolution = con_res[j]) #5
      u <- UMAPPlot(object =sub , label=T) + 
        NoLegend() + 
        labs(title=paste("Dimention = ", con_dim[i], ",  Resolution = ", con_res[j]), sep="") + 
        theme(plot.title = element_text(hjust = 0.5))
      plot(u) 
    }
  }
  dev.off()
}





###### dimensionとresolutionのパラメータを変えながら，クラスタリングを行ってみる関数
fun.cluster.condition <- function(obj, filename){
  DefaultAssay(obj) <- "integrated"
  
  obj <- RunPCA(obj, verbose = FALSE)
  
  pdfname <- paste("ClusterFinding-", filename, ".pdf", sep="")
  pdf(pdfname, width = 6, height = 6)
  
  for (i in 1:length(con_dim)){
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:con_dim[i])
    obj <- FindNeighbors(object = obj, dims = 1:con_dim[i])
    for (j in 1:length(con_res)){
      obj <- FindClusters(object = obj, resolution = con_res[j])
      u <- UMAPPlot(object = obj , label=T) + 
        NoLegend() + 
        labs(title=paste("Dimention = ", con_dim[i], ",  Resolution = ", con_res[j]), sep="") + 
        theme(plot.title = element_text(hjust = 0.5))
      plot(u) 
    }
  }
  dev.off()
}

# Execution
con_dim = c(30,35) #綺麗でない
con_res = seq(0.1, 0.4, by=0.05)
fun.cluster.condition(wh.woRBC.sct.v3, "WHwoRBC_SCTv3(CCscore)") 
# candidates -> 
con_dim = c(27,28,29,30) #31,36はだめ
con_res = seq(0.05, 0.4, by=0.05)
fun.cluster.condition(wh.woRBC.sct.v3, "WHwoRBC_SCTv3(CCscore)") 