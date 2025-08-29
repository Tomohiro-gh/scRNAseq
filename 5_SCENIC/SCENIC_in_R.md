## SCENIC in R
https://github.com/aertslab/SCENIC

### part1) Introduction & Setup
https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html

### part2) Running SCENIC : 
https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html


database [here](https://resources.aertslab.org/cistarget/databases/old/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/)


- mm9-500bp-upstream-7species.mc9nr.feather
- mm9-tss-centered-10kb-7species.mc9nr.feather

## SENIC set up in R (24/10/24)

Tutorial codes are below:
```r
library(SCENIC)
org <- "mgi" # or hgnc, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Mouse brain" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
```
1) Error1: 
Trouble shooting: 


2) Error2: motifAnnotations_mgi not found
`initializescenic` でcisTargets databeseのファイルを置き換える :  https://github.com/aertslab/SCENIC/issues/364

->
```r
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
```

```r
org <- "mgi" # or hgnc, or dmel
dbDir <- "cisTarget_databases/mm10"
myDatasetTitle <- "SCENIC" # choose a name for your analysis
data(defaultDbNames)
defaultDbNames$mgi[1] <- "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
defaultDbNames$mgi[2] <- "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather"
defaultDbNames
data(list="motifAnnotations_mgi_v9", package="RcisTarget")
motifAnnotations_mgi <- motifAnnotations_mgi_v9
scenicOptions <- initializeScenic(org = org,
                                  dbDir = dbDir,
                                  dbs = defaultDbNames[["mgi"]],
                                  datasetTitle = myDatasetTitle,
                                  nCores = cores)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$seed <- 123
```


