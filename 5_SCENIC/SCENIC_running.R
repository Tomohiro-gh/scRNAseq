

# https://rdrr.io/github/aertslab/SCENIC/f/vignettes/SCENIC_Running.Rmd

## Initialize SCENIC settings

##  The parameters that need to be specified in all runs is the organism (mgi for mouse, hgnc for human, or dmel for fly), and the directory where the RcisTarget databases are stored (you may create a link in the current directory to avoid duplicating them, e.g. in linux: system("ln -s ~/path/to/dbs databases")).

library(SCENIC)
org <- "mgi" # or hgnc, or dmel
# dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Mouse brain" # choose a name for your analysis

dbDir = "/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/SCENIC/RcisTarget_database/Mouse/mm9_r45_mc9nr" ##mm9のふるいやつ

## https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/
## ここらかダウンロードした．新しいdbではエラーが出る．

## Error in initializeScenic(org = org, dbDir = dbDir, dbs = dbs, datasetTitle = myDatasetTitle,  : 
# The following RcisTarget databases were not found: 
#   - /Users/tomohiro/Dropbox/FukuharaLab_Res/Database/SCENIC/RcisTarget_database/Mouse/mm10_mc_v10_clust//mm9-500bp-upstream-7species.mc9nr.feather 
# - /Users/tomohiro/Dropbox/FukuharaLab_Res/Database/SCENIC/RcisTarget_database/Mouse/mm10_mc_v10_clust//mm9-tss-centered-10kb-7species.mc9nr.feather
# Make sure the arguments 'dbDir' and 'dbs' are correct.
# 

data(defaultDbNames)

dbs <- defaultDbNames[[org]]

##  重要：https://github.com/aertslab/SCENIC/issues/364
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9

scenicOptions <- 
  initializeScenic(
    org=org, 
    dbDir=d
    dbs=dbs, 
    datasetTitle=myDatasetTitle,
    nCores=10) 

scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123


# Modify if needed
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
# Databases:
# scenicOptions@settings$dbs <- c("mm9-5kb-mc8nr"="mm9-tss-centered-5kb-10species.mc8nr.feather")
# scenicOptions@settings$db_mcVersion <- "v8"

# Save to use at a later time...
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
