# Single Cell network analysis

## SCENIC (Single Cell rEgulatory Network Inference and Clustering) (for scRNAseq. data)
#### R or Python is avairable 
to install in R

https://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Setup.html

```r
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

```

さらに，databaseの構築が必要

今回はこちらを使用

Mouse：https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/

Human：https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/

```sh
# Terminalで．
cd xxx

For mouse
feather_database_url='https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
wget "${feather_database_url}"
feather_database_url2='https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather'
wget "${feather_database_url2}"
```


Workflow 

https://rdrr.io/github/aertslab/SCENIC/f/vignettes/SCENIC_Running.Rmd


------------------
## scGRN (for multi-omics)
web based 

https://bio.liclab.net/scGRN/

------------------
## CellOracle (for scRNAseq. data)
Documentation: https://morris-lab.github.io/CellOracle.documentation/index.html

Github: https://github.com/morris-lab/CellOracle

Article: https://www.nature.com/articles/s41586-022-05688-9


------------------
## DeepMAPS (for multi-omics)
Docker:  https://github.com/OSU-BMBL/deepmaps

Webserver: https://bmblx.bmi.osumc.edu/

Article: https://www.nature.com/articles/s41467-023-36559-0


------------------
## SCENIC+ (for multi-omics)
https://scenicplus.readthedocs.io/en/latest/install.html

Article: https://www.nature.com/articles/s41592-023-01938-4

To install To install SCENIC+ run in conda environment
```sh
cd Documents
conda create --name scenicplus python=3.8
conda activate scenicplus
git clone https://github.com/aertslab/scenicplus
cd scenicplus
pip install -e .
```


Tutorials: https://scenicplus.readthedocs.io/en/latest/tutorials.html#
