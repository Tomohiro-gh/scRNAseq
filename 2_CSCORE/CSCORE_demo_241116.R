library(CSCORE) #install_github("ChangSuBiostats/CS-CORE")
library(WGCNA) #BiocManager::install("WGCNA")
library(impute)
library(Seurat)

#> https://changsubiostats.github.io/CS-CORE/index.html
#> https://changsubiostats.github.io/CS-CORE/articles/CSCORE.html
#> 
#> 
obj <- readRDS("/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp155_reanalysis_of_2022DC_LCangiogenesis/3_ClusterAnnotation_Dim26Res0.5/26782Cells_Dim26Res0.5_annotated.rds")

# Use the original UMI counts stored in Assay 'RNA'
DefaultAssay(pbmc) <- 'RNA'

pc = obj[ , obj$CellAnnotation_Fine %in% 'Pericyte']
  pc 
  # An object of class Seurat 
  # 20478 features across 610 samples within 2 assays 
  # Active assay: RNA (18478 features, 0 variable features)
  # 3 layers present: counts, data, scale.data
  # 1 other assay present: integrated
  # 3 dimensional reductions calculated: pca, umap, tsne

mean_exp = rowMeans(pc@assays$RNA@counts/pc$nCount_RNA)
genes_selected = names(sort.int(mean_exp, decreasing = T))[1:5000]

pc_wo <- pc[, pc$State == "WO"]


CSCORE_result <- CSCORE(pc_wo, genes = genes_selected)


#> 4. Downstream analysis on the co-expression network
#> 4.1 Extract co-expressed gene module
#> 
# Obtain CS-CORE co-expression estimates
CSCORE_coexp <- CSCORE_result$est

# Obtain BH-adjusted p values
CSCORE_p <- CSCORE_result$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp[p_matrix_BH > 0.05] <- 0


# Compute the adjacency matrix based on the co-expression matrix
adj = WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 1)
# Compute the topological overlap matrix
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
# Run hierarchical clustering as in the WGCNA workflow
hclust_dist = hclust(as.dist(dissTOM), method = "average") 
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                     distM = dissTOM, 
                                     deepSplit = 2,
                                     pamRespectsDendro = FALSE,
                                     minClusterSize = 10)
# For more instructions on how to tune the parameters in the WGCNA workflow,
# please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

names(memb) = genes_selected
memb_tab <- table(memb)
module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))

module_list[[1]]

## > ここからoriginal
names(module_list) <- paste("GeneModule",
                            seq(1, length(module_list), by=1),
                            sep = "_")
  #> do not run
  names(module_list)[1]
  module_list[module_list=="Pde5a"]

  
## ここから特定の遺伝子を探してみる
#> do not run
for(i in 1:length(module_list)){
    if(any(module_list[[i]] == "Tomm7") == TRUE){
      print(names(module_list[i]))
    }
  }


FUN.search.genes.oi <- function(list, gene_oi){
  for(i in 1:length(module_list)){
    if(any(module_list[[i]] == gene_oi) == TRUE){
      print(names(module_list[i]))
    }
  }
}
  FUN.search.genes.oi(module_list, "Gsk3a")
  
  
  FUN.search.genes.oi(module_list, "Pdgfrb")
  FUN.search.genes.oi(module_list, "Col6a1")
  
saveRDS(module_list, "module_list_PC_WO.rds")
  