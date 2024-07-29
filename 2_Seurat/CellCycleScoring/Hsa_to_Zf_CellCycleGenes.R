# CellCycleGeneを抜き出す

library(Seurat) # cc.genesがロードされる

dir = "/Users/tomohiro/Dropbox/FukuharaLab_Res/Zebrafish/Database_zf/ensembl_biomart/BiomaRt_ensembl_zf_hsa_ortholog_list.xlsx"

ZftoHsa <- read.xlsx(dir, rowNames = FALSE, startRow = 1)　
head(ZftoHsa)
# ensembl_gene_id external_gene_name hsapiens_homolog_ensembl_gene hsapiens_homolog_associated_gene_name

cc.genes
hsa.s.genes <- cc.genes$s.genes
  length(hsa.s.genes) #43
hsa.g2m.genes <- cc.genes$g2m.genes
  length(hsa.g2m.genes) #54

zf.s.genes <- ZftoHsa %>%
  filter(hsapiens_homolog_associated_gene_name %in% hsa.s.genes)  %>%
    pull(external_gene_name)
length(zf.s.genes) #42

zf.g2m.genes <- ZftoHsa %>%
  filter(hsapiens_homolog_associated_gene_name %in% hsa.g2m.genes)  %>%
    pull(external_gene_name)
length(zf.g2m.genes) #53

# listに名前をつけて保存
cc.genes.zf <- setNames(list(zf.s.genes, zf.g2m.genes), 
                        eval(substitute(alist("zf.s.genes", "zf.g2m.genes"))))

saveRDS(cc.genes.zf, file="CellCycleGenes_ZF.rds")



#saveして移動
cc.genes.zf <- readRDS("/Users/tomohiro/Dropbox/FukuharaLab_Res/Zebrafish/Database_zf/CellCycleGenes_ZF.rds")
zf.s.genes <- cc.genes.zf$zf.s.genes
zf.g2m.genes <- cc.genes.zf$zf.g2m.genes