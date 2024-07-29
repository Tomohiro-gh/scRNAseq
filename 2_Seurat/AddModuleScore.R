## AddModuleScore Function
library(Seurat)


## Official
AddModuleScore(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = "Cluster",
  seed = 1,
  search = FALSE,
  ...
)


## Usage
object
#Seurat object

features
#A list of vectors of features for expression programs; each entry should be a vector of feature names

pool
#List of features to check expression levels against, defaults to rownames(x = object)

nbin
#Number of bins of aggregate expression levels for all analyzed features

ctrl
#Number of control features selected from the same bin per analyzed feature

k
#Use feature clusters returned from DoKMeans

assay
#Name of assay to use

name
#Name for the expression programs; will append a number to the end for each entry in features (eg. if features has three programs, the results will be stored as name1, name2, name3, respectively)

seed
#Set a random seed. If NULL, seed is not set.

search
#Search for symbol synonyms for features in features that don't match features in object? Searches the HGNC's gene names database; see UpdateSymbolList for more details

...
#Extra parameters passed to UpdateSymbolList




#### 
if (FALSE) {
  data("pbmc_small")
  cd_features <- list(c(
    'CD79B',
    'CD79A',
    'CD19',
    'CD180',
    'CD200',
    'CD3D',
    'CD2',
    'CD3E',
    'CD7',
    'CD8A',
    'CD14',
    'CD1C',
    'CD68',
    'CD9',
    'CD247'
  ))
  pbmc_small <- AddModuleScore(
    object = pbmc_small,
    features = cd_features,
    ctrl = 5,
    name = 'CD_Features'
  )
  

