### set up at 05/16/23 

## liana ++ 
## https://saezlab.github.io/liana/articles/liana_devel.html#install-liana

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

library(devtools)

Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE) # ignore warning from iTALK
BiocManager::install("ComplexHeatmap") # required for Connectome
devtools::install_github('saezlab/OmnipathR@ff3ad88e3915747e1b557bf44ac5396f9525dd7e') # install 4.0 version of OmnipathR

# install tools
devtools::install_github("sqjin/CellChat")
devtools::install_github('msraredon/Connectome', ref = 'master')
devtools::install_github("Coolgenome/iTALK", build_vignettes = FALSE)
# A modified version of SingleCellSignalR (SCA) that enables external resources
devtools::install_github(repo = "saezlab/SingleCellSignalR_v1",
                         subdir = "SingleCellSignalR")

# Finally, install LIANA
devtools::install_github('saezlab/liana')


### ここはコマンドラインで
## python
## conda create -n liana_env
## conda activate liana_env
## conda install -c anaconda python=3.8.5
## pip install squidpy

# liana.ymlファイルをダンロード -> home directoryへ．
# https://github.com/saezlab/liana/blob/master/liana_env.yml

## home directoryで，
## conda env create -f liana_env.yml

#3 Termianalで
# path.package("liana")
## "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/liana"
# cd "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/liana/"
# git clone https://github.com/saezlab/NATMI


# If needed, use the following to locate the liana package
system.file(package = "liana")


## 最終的に，
## condaで
## liana_env -> conda createで作ったもの
## liana_env2 -> liana_env.ymlを読み込ませたもの
