## 24/03/15

##SCENIC

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::version() #[1] ‘3.18’
# If your bioconductor version is previous to 4.0, see the section bellow

## Required
BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
BiocManager::install(c("RcisTarget"), version = 3.17)

library(AUCell)
# library(RcisTarget) #biocmanager 3.18’だとインストール不可
library(GENIE3)

#https://github.com/aertslab/RcisTarget/issues/8
library(remotes)
remotes::install_github("aertslab/RcisTarget")
library(RcisTarget)


## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
## rbokehもversion3.18じゃだめ
remotes::install_github("bokeh/rbokeh")

# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)


if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
packageVersion("SCENIC") #[1] ‘1.3.1’


dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather")


for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}
## OK



