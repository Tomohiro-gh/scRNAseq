## Liana tutorial
## https://saezlab.github.io/liana/articles/liana_tutorial.html

library(tidyverse)
library(magrittr)
library(liana)


#### CCC Resources
show_resources()
# [1] "Default"    "Consensus" "Baccin2019"  "CellCall"   "CellChatDB" "Cellinker"       
# [7] "CellPhoneDB" "CellTalkDB"  "connectomeDB2020" "EMBRACE" "Guide2Pharma"  "HPMR"            
# [13] "ICELLNET"  "iTALK" "Kirouac2010"  "LRdb" "Ramilowski2015"  "OmniPath"        
# [19] "MouseConsensus"  

#### CCC Methods
show_methods()
# [1] "connectome" "logfc" "natmi"  "sca"  "cellphonedb"  "cytotalk" "call_squidpy"   
# [8] "call_cellchat"   "call_connectome" "call_sca"   "call_italk"   "call_natmi"     


liana_path <- system.file(package = "liana")
testdata <-
  readRDS(file.path(liana_path , "testdata", "input", "testdata.rds"))

testdata %>% dplyr::glimpse()

# Run liana
liana_test <- liana_wrap(testdata)
#> Warning in exec(output, ...): 3465 genes and/or 0 cells were removed as they had
#> no counts!
#> Warning: `invoke()` is deprecated as of rlang 0.4.0.
#> Please use `exec()` or `inject()` instead.
#> This warning is displayed once per session.

# Liana returns a list of results, each element of which corresponds to a method
liana_test %>% dplyr::glimpse()




##### Aggregate and Obiain Consensus Ranks
# We can aggregate these results into a tibble with consensus ranks
liana_test <- liana_test %>%
  liana_aggregate()

dplyr::glimpse(liana_test)


## Simple DotPlot
liana_test %>%
  liana_dotplot(source_groups = c("B"),
                target_groups = c("NK", "CD8 T", "B"),
                ntop = 20)


## Frequency Heatmap
liana_trunc <- liana_test %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)


## Frequency Chord diagram
if(!require("circlize")){
  install.packages("circlize", quiet = TRUE,
                   repos = "http://cran.us.r-project.org")
}
p <- chord_freq(liana_trunc,
                source_groups = c("CD8 T", "NK", "B"),
                target_groups = c("CD8 T", "NK", "B"))





##### Run any method of choice.
# Load Sce testdata
sce <- readRDS(file.path(liana_path , "testdata", "input", "testsce.rds"))

# RUN CPDB alone
cpdb_test <- liana_wrap(sce,
                        method = 'cellphonedb',
                        resource = c('CellPhoneDB'),
                        permutation.params = list(nperms=100,
                                                  parallelize=FALSE,
                                                  workers=4),
                        expr_prop=0.05)

# identify interactions of interest
cpdb_int <- cpdb_test %>%
  # only keep interactions with p-val <= 0.05
  filter(pvalue <= 0.05) %>% # this reflects interactions `specificity`
  # then rank according to `magnitude` (lr_mean in this case)
  rank_method(method_name = "cellphonedb",
              mode = "magnitude") %>%
  # keep top 20 interactions (regardless of cell type)
  distinct_at(c("ligand.complex", "receptor.complex")) %>%
  head(20)



# Plot toy results
cpdb_test %>%
  # keep only the interactions of interest
  inner_join(cpdb_int, 
             by = c("ligand.complex", "receptor.complex")) %>%
  # invert size (low p-value/high specificity = larger dot size)
  # + add a small value to avoid Infinity for 0s
  mutate(pvalue = -log10(pvalue + 1e-10)) %>% 
  liana_dotplot(source_groups = c("c"),
                target_groups = c("c", "a", "b"),
                specificity = "pvalue",
                magnitude = "lr.mean",
                show_complex = TRUE,
                size.label = "-log10(p-value)")



# Run liana re-implementations with the CellPhoneDB resource
complex_test <- liana_wrap(testdata,
                           method = c('natmi', 'sca', 'logfc'),
                           resource = c('CellPhoneDB'))

complex_test %>% liana_aggregate()


##### Call liana with overwritten default settings
# define geometric mean
geometric_mean <- function(vec){exp(mean(log(vec)))}

# Overwrite default parameters by providing a list of parameters
liana_test <- liana_wrap(testdata,
                         method = c('cellphonedb', 'sca'),
                         resource = 'Consensus',
                         permutation.params = 
                           list(
                             nperms = 10 # here we run cpdb it only with 10 permutations
                           ),
                         complex_policy="geometric_mean"
)

# This returns a list of results for each method
liana_test %>% dplyr::glimpse()