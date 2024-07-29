library(Seurat)

## Cell cycle scoring
## original code: 220524_1_integration_standard_allcells.R

##################  CellCycle Scoreを追加  ####################################
# code from: https://github.com/satijalab/seurat/issues/3692

# importing cellcycle gene list in zf
cc.genes.zf <- readRDS("/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/CellCycleGenes_ZF.rds")
cc.genes.zf <- setNames(object=cc.genes.zf, c("zf.s.genes", "zf.g2m.genes"))
zf.s.genes <- cc.genes.zf$zf.s.genes
zf.g2m.genes <- cc.genes.zf$zf.g2m.genes


cc.genes.zf 


