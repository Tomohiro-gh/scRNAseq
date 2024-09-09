library(Seurat)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
######### Dim の決定 : elbow plot
fun.elb.plot <- function(object,filename){
  elb <- ElbowPlot(object, ndims = 50) + 
    labs(title=paste0("ElbowPlot : ", filename)) + 
    theme(plot.title = element_text(hjust = 0.5))
  plot(elb)
  
  ggsave(
    paste0("ClusteringCondition_Elbow_",filename,".png"),
    width=7, height=5, dpi=300)
}
# Execution
# fun.elb.plot(woRBC8898_int, "1355Cells")
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>> Dim Heatmap >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 2個ずつ書いていく
fun.Dim.Heat <- 
  function(SeuratObject, dimention, Title){
    pdfname <- 
      paste0("ClusteringCondition_DimHeat_", Title, ".pdf")
    
  pdf(pdfname, width = 10, height = 8)
  DimHeatmap(
    SeuratObject, 
    dims = 1:20, 
    cells = 500, 
    balanced = TRUE)
  mtext(side = 3, line=0.2, outer=T, text = Title, cex=1)
  par(oma = c(0, 0, 1.5, 0)) 
  options(repr.plot.width = 15, repr.plot.height = 15)
    #Dimを21以上に指定した場
  if(dimention >= 21){
    DimHeatmap(
      SeuratObject, 
      dims = 21:dimention, 
      cells = 500, 
      balanced = TRUE)
    mtext(side = 3, line=0.2, outer=T, text = Title, cex=1)
    par(oma = c(0, 0, 1.5, 0)) 
    options(repr.plot.width = 15, repr.plot.height = 15)
    dev.off()
    
  }else{dev.off()}
}
# example
# example: fun.Dim.Heat(ECPC, 1:20, "ECPC, standard integration")
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>> Cluster finding >>>>>>>>>>>>>>>>>>>>> 05/23
# This function is requiered for predefined dim, res range like above
# PCAまで行っているものが対象
# Usage: fun.cluster.finding.UMAP(sub, filename)
# Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/SeuratFunctions.R
#con_dim = c(26,31,38,46) # variable depending on condition
# con_res = seq(0.1, 0.3, by=0.05) # variable depending on condition
fun.cluster.finding.UMAP <- function(sub, filename){
  pdfname <- paste0("ClusteringConditionUMAP_", filename, ".pdf")
  
  pdf(pdfname, width = 6, height = 6)
  
  for (i in 1:length(con_dim)){
    sub <- RunUMAP(sub, dims = 1:con_dim[i]) #3
    sub <- FindNeighbors(object = sub, dims = 1:con_dim[i]) #4
    for (j in 1:length(con_res)){
      sub <- FindClusters(object = sub, resolution = con_res[j]) #5
      u <- UMAPPlot(object =sub , label=T) + 
        NoLegend() + 
        labs(title=paste0("Dimention = ", con_dim[i], ",  Resolution = ", con_res[j])) + 
        theme(plot.title = element_text(hjust = 0.5))
      plot(u) 
    }
  }
  dev.off()
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
