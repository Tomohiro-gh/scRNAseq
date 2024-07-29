### suerat objectをdataframeへ変換する


# Load data.frame obj
pbmc <- readRDS("data/pbmc_2k_v3_df.rds")
identity <- readRDS("data/pbmc_2k_v3_Seurat_Idents.rds")

features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

# Subset data.frame
pbmc <- pbmc[,features]

# Add cell ID and identity classes
pbmc$Cell <- rownames(pbmc)
pbmc$Idents <- identity

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Idents"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")
head(pbmc, 10)