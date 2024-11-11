library(SoupX)

sc = load10X('/Volumes/Ishii@FukuharaLab-2/Exp190_agedHeartEC_scRNAseq/1_cellranger_count_v1/190-1_3month-1/outs')
sc = autoEstCont(sc)
out = adjustCounts(sc)

