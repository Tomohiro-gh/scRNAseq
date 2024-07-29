library(Seurat)

## Cell cycle scoring
## original code: 220524_1_integration_standard_allcells.R

##################  CellCycle Scoreを追加  ####################################
# code from: https://github.com/satijalab/seurat/issues/3692

# importing cellcycle gene list in zf
cc.genes.zf <- readRDS("/Users/tomohiro/Dropbox/FukuharaLab_Res/Database/CellCycleGenes_ZF.rds")
cc.genes.zf <- readRDS("/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/CellCycleScoring/CellCycleGenes_ZF.rds")


cc.genes.zf <- setNames(object=cc.genes.zf, c("zf.s.genes", "zf.g2m.genes"))
zf.s.genes <- cc.genes.zf$zf.s.genes
zf.g2m.genes <- cc.genes.zf$zf.g2m.genes

saveRDS(cc.genes.zf, "/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/CellCycleScoring/CellCycleGenes_ZF.rds")


##CellCycle scoring ###### * lapplyはダメ．metadataの確認ができない．
for (i in 1:length(woRBC8898.list)) {
  woRBC8898.list[[i]] <- NormalizeData(woRBC8898.list[[i]], verbose = FALSE)
  woRBC8898.list[[i]] <- CellCycleScoring(woRBC8898.list[[i]],
                                          s.features = zf.s.genes,
                                          g2m.features = zf.s.genes,
                                          set.ident = TRUE)
  woRBC8898.list[[i]]$cc_difference <- woRBC8898.list[[i]]$S.Score - woRBC8898.list[[i]]$G2M.Score
}



