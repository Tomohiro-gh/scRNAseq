#> Seurat install 
install.packages("remotes")
library(remotes)
#> Install previous versions of Seurat
#> https://satijalab.org/seurat/articles/install_v5.html#install-previous-versions-of-seurat
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))



