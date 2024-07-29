## Seuratの環境構築編

## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Functionをlocal環境へload
# Load a script only to my environment
Dir="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/LoadingFunctions/"
scripts <- 
  c("function_Graphs_v2.R",
    "function_FindAllMarkers.R", 
    "function_ClusteringCondition.R",
    "function_check_Clustering_UMAPplot.R")

.myfunc.env = new.env()
for(sc in scripts){
  sys.source(paste0(Dir, sc), envir = .myfunc.env ) # localにload
}

sys.source("/Users/tomohiro/Dropbox/FukuharaLab_Res/Experiment/Exp174_tubalmuris_MCA_EC/Heart/Heart_MCA_for_Umapplot.R", envir = .myfunc.env)

attach(.myfunc.env)
## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<