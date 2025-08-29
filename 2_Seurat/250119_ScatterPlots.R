library(tidyplots)
library(matrix)


obj <- readRDS('/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp190_AgedHeratEC_scRNAseq/Part4_EC_3,12,21month/P4-0_ECintegration_3,12,21month/HeartEC_3,12,21month_v1.rds')
obj


## scatter plot
GetAssayData(obj, )
counts <- obj[["RNA"]]$data
head(counts)
t(counts) %>% head

#> sample
cnt <- obj[["RNA"]]$data
df_cnt <- t(as.matrix(cnt[c("Tgfb2", "Smad7"),]))
df_mt <- obj@meta.data |> dplyr::select(Age)
df_tmp <- cbind(df_mt, df_cnt)
  head(df_tmp)



#> ２つの遺伝子を抽出し，dataframeにする
#> 
FUN.Scatter.2Genes <- function(seurat_obj, GeneX, GeneY, Ident_oi){
  countlayer <- seurat_obj[["RNA"]]$data
  df_count <- t(as.matrix(countlayer[c(GeneX, GeneY),]))
  df_meta <- seurat_obj@meta.data |> dplyr::select(as.name(Ident_oi))
  tmp <- cbind(df_meta, df_count)
    print(head(tmp))
  g1 <- ggplot(data = tmp,
           aes_(x=as.name(GeneX), y=as.name(GeneY), color = as.name(Ident_oi)))+
    geom_point()
  plot(g1)
}
#> Execution
FUN.Scatter.2Genes(obj, "Pecam1", "Kdr", Ident_oi = "Age")


df_tmp |>
  tidyplot(x=Tgfb2, y=Smad7, color = Age) |>
  add_data_points(size = 2, alpha = 0.7) 
plot(t1)
df_tmp |>
  tidyplot(x=Tgfb2, y=Smad7) |>
  add_data_points() 
plot(t1)


#> tidyで書いてみる
FUN.Scatter.2Genes.tidy <- function(seurat_obj, GeneX, GeneY, Ident_oi){
  countlayer <- seurat_obj[["RNA"]]$data
  df_count <- t(as.matrix(countlayer[c(GeneX, GeneY), ]))
  df_meta <- seurat_obj@meta.data |> dplyr::select(as.name(Ident_oi))
  tmp <- cbind(df_meta, df_count)
  colnames(tmp) <- c("group", "geneX", "geneY")
    print(head(tmp))
  t1 <- tmp %>% 
    tidyplot(x=geneX, y=geneY, color = group) %>% 
    add_data_points() 
  plot(t1)
}
#> Execution
FUN.Scatter.2Genes.tidy(obj, "Pecam1", "Kdr", Ident_oi = "Age")
## Ave Exp vs Time
A = "Age"
as.name(A)


