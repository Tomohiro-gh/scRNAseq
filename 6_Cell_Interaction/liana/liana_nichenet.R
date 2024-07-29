library(tidyverse)
library(liana)
library(nichenetr)
library(Seurat)
library(ggrepel)
library(cowplot)

hnscc_expression = readRDS(
  "/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/hnscc_expression.rds")

expression = hnscc_expression$expression
expression[1:3,1:3]

sample_info = hnscc_expression$sample_info