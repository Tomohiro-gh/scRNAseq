#> cluster profiler

library(clusterProfiler)
library(re)


## Gene id conversion 
## 色々あるが，clusterProfilerのbitr functionを使えばいい
## https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html#reactomepa-supported-organisms

gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)


#>GO
#>https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

#>GSEA
#>
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

#>KEGG
#>https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
#> 使えるorganismの検索
search_kegg_organism('ece', by='kegg_code')

kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

ekg <- enrichKEGG(gene = names(Named_logFC_list), 
                  universe = all.genes.df$entrezgene_id,
                  organism = "mmu",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  minGSSize = 5)

## 7.5 KEGG module gene set enrichment analysis
## https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-kegg.html
mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa',
                 pvalueCutoff = 1)


#>Reactome
#> https://yulab-smu.top/biomedical-knowledge-mining-book/reactomepa.html
library(ReactomePA)
data(geneList, package="DOSE"
)
de <- names(geneList)[abs(geneList) > 1.5]
x <- enrichPathway(gene=de, pvalueCutoff = 0.05, readable=TRUE)

y <- gsePathway(geneList, 
                pvalueCutoff = 0.2,
                pAdjustMethod = "BH", 
                verbose = FALSE)
