## https://rdrr.io/github/aertslab/SCENIC/f/vignettes/SCENIC_Running.Rmd

library(SCENIC)
library(SCopeLoomR)


wd = "/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/SCENIC/demo_Brain"
setwd(wd)


loomPath <- system.file(package="SCENIC", "examples/mouseBrain_toy.loom")

## Open the loom file and load the expression matrix (and cell annotation if available)
loom <- open_loom(loomPath)
exprMat <- get_dgem(loom)
cellInfo <- get_cell_annotation(loom)
close_loom(loom)
  dim(exprMat)


## Cell info/phenodata
# cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)
cellInfo <- data.frame(cellInfo)
cbind(cellInfo, table(cellInfo$CellType))
dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")
# Color to assign to the variables (same format as for NMF::aheatmap)
colVars <- list(CellType=c("microglia"="forestgreen", 
                           "endothelial-mural"="darkorange", 
                           "astrocytes_ependymal"="magenta4", 
                           "oligodendrocytes"="hotpink", 
                           "interneurons"="red3", 
                           "pyramidal CA1"="skyblue", 
                           "pyramidal SS"="darkblue"))
colVars$CellType <- 
  colVars$CellType[
    intersect(names(colVars$CellType),
              cellInfo$CellType)]

saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))


## Initialize SCENIC settings

library(SCENIC)
org <- "mgi" # or hgnc, or dmel
## (mgi for mouse, hgnc for human, or dmel for fly),

dbDir <- "/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/SCENIC/RcisTarget_database/Mouse/mm9_r45_mc9nr" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Mouse brain" # choose a name for your analysis

##  重要：https://github.com/aertslab/SCENIC/issues/364
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

data(defaultDbNames)
dbs <- defaultDbNames[[org]]

scenicOptions <- 
  initializeScenic(
    org = org, 
    dbDir= dbDir, 
    dbs = dbs,
    datasetTitle = myDatasetTitle,
    nCores = 10) 

# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## 呼び出す時:
 


## Co-expression network


## Gene filter/selection
#### To run GENIE3/GRNBoost we recommend to apply soft gene filter, to remove genes that are expressed either at very low levels or in too few cells.
#### Filter by the total number of reads per gene.
genesKept <- 
  geneFiltering(
    exprMat,
    scenicOptions=scenicOptions,
    minCountsPerGene=3*.01*ncol(exprMat),
    minSamples=ncol(exprMat)*.01)
  length(genesKept)

## check 
interestingGenes <- c("Sox9", "Sox10", "Dlx5")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)] #取り逃がしたGOIがなければOK

#### filter the expression matrix to contain only these
exprMat_filtered <- exprMat[genesKept, ]
  dim(exprMat_filtered)
  rm(exprMat)
  

# Correlation
runCorrelation(exprMat_filtered, scenicOptions)
  

# GENIE3
#### Since GENIE3 is based on a Random Forest approach, each time it is run the results will be slightly different. The higher the number of trees used (ntrees), the lower the variability. We recommend to use set.seed to reproduce exact results in multiple runs. For more details, check ?GENIE3 (GENIE3 help) or ?runGenie3 (SCENIC wrapper for GENIE3).


# Optional: add log (if it is not logged/normalized already)
exprMat_filtered <- log2(exprMat_filtered+1) 

# Run GENIE3
runGenie3(exprMat_filtered, scenicOptions)
  

# Build and score the GRN (runSCENIC_…)

#### Re-load the expression matrix if necessary:
   loom <- open_loom(loomPath)
   exprMat <- get_dgem(loom)
   close_loom(loom)
  # # Optional: log expression (for TF expression plot, it does not affect any other calculation)
   exprMat_log <- log2(exprMat+1)
   dim(exprMat)

#### Run the remaining steps using the wrapper functions:
library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

# For a very quick run: 
# coexMethod=c("top5perTarget")
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
# save...

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)

saveRDS(scenicOptions, file="int/scenicOptions.Rds") # To save status

