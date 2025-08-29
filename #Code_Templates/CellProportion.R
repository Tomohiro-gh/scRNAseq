##> Cellproportion plots
##> 250117_cellproportion_HeartECs.R をreferenceに
library(Seurat)
library(ggplot2)
library(dplyr)
library(openxlsx)
library(scales)

wd = '/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp190_AgedHeratEC_scRNAseq/Part4_EC_3,12,21month/P4-1_CellProportion'
setwd(wd)

## ------------------------------------------------------------------
## Load a script only to my environment
.myfunc.env = new.env()
## for scRNAseq
Dir="/Users/tomohiro/Library/CloudStorage/Dropbox/GitHub/scRNAseq/#LoadingFunctions/"
(scripts = list.files(Dir,pattern = ".R$"))
for(sc in scripts){
  sys.source(paste0(Dir, sc), envir = .myfunc.env ) # localにload
  }
attach(.myfunc.env)
## -----------------------------------------------------------------

#> Heart EC 3,12,21 month
obj <- readRDS('/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp190_AgedHeratEC_scRNAseq/Part4_EC_3,12,21month/P4-0_ECintegration_3,12,21month/HeartEC_3,12,21month_v1.rds')
  obj
  # An object of class Seurat 
  # 32285 features across 50732 samples within 1 assay 
  # Active assay: RNA (32285 features, 2000 variable features)

  colnames(obj@meta.data) %>% data.frame()
  # 1            orig.ident
  # 2            nCount_RNA
  # 3          nFeature_RNA
  # 4            percent.mt
  # 5            SampleName
  # 6                Tissue
  # 7                   Age
  # 8            percent.rb
  # 9             pANN_0.25
  # 10     DF_classfication
  # 11       EC_int_cca_sub
  # 12 EC_int_cca_sub_x_Age
  # 13      RNA_snn_res.0.2
  # 14      RNA_snn_res.0.4
  # 15      RNA_snn_res.0.6
  # 16      RNA_snn_res.0.8
  # 17        RNA_snn_res.1
  # 18      RNA_snn_res.1.2
  # 19      RNA_snn_res.1.4
  # 20      RNA_snn_res.1.6
  # 21      RNA_snn_res.1.8
  # 22        RNA_snn_res.2
  # 23      seurat_clusters
  # 24     umap.cca.cluster
  # 25              ECannot
  # 26         ECannot_Fine
  
  
#> function ----------------------------------------------
#> colnameは　Age, Ratio, Celltypeで固定
FUN.ProportionPlot <- function(DataFrame, ClusterColors, FileName, GraphLabel=NULL){
  
  #> themes
  theme_bar <- theme(plot.title = element_text(hjust = 0.5, size = 16), 
                     plot.subtitle = element_text(size = 14, hjust = 0.5),
                     axis.title.x = element_text(hjust= 0.5, size = 16, face = "bold"),
                     axis.title.y = element_text(size = 16, face = "bold"),
                     axis.text.x = element_text(colour= "black", size= 14,  face ="bold"),
                     axis.text.y = element_text(colour = "black", size= 14),
                     axis.line=element_line(colour = "black", linewidth=1.0), # 軸の色
                     axis.ticks=element_line(colour = "black", linewidth = 1.0), #軸の線
                     axis.ticks.length = unit(3, "mm"),
                     panel.grid.major = element_blank(), # no main grid line
                     panel.grid.minor = element_blank(), # no sub grid line
                     panel.background = element_blank(), # white background
                     text = element_text(size = 16, face = "bold"))
  
  #> 1 横軸： stage, 縦軸ratio - 数値を入れないversion
  g1 <- ggplot(DataFrame, aes(x = Age, y = Ratio, fill = Celltype)) + 
    geom_bar(width = 0.6, stat = "identity", position = "fill") + 
    scale_y_continuous(expand=c(0,0), labels = percent) +
    scale_fill_manual(values = ClusterColors) +
    labs(x="", y="",title = "") + 
    scale_fill_discrete(name=GraphLabel) +
    theme_bar
      plot(g1)
        ggsave(paste0(FileName, "_1_Proportion_tate.png"),
               width= 8, height = 8, dpi=300)
  ## 2 横軸： stage, 縦軸ratio - ラベルあり
  g2 <- g1 + 
     geom_text(data = DataFrame,
            aes(y = Ratio, label = Percent_show),
            position = position_fill(vjust = 0.5))
      plot(g2)
        ggsave(paste0(FileName, "_2_Proportion_tate_label.png"),
             width= 8, height = 8, dpi=300)
  
  ## 3 横向きグラフ
  g3 <- ggplot(DataFrame, aes(x = Ratio, y = Age, fill = Celltype)) + 
    geom_bar(width = 0.6, stat = "identity", position = "fill") + 
    scale_x_continuous(expand=c(0,0), labels = percent) +
    scale_fill_manual(values = ClusterColors) +
    labs(x="", y="",title = "") + 
    theme_bar + 
    scale_fill_discrete(name=GraphLabel) 
      plot(g3)
        ggsave(paste0(FileName, "_3_Proportion_yoko.png"),
             width= 12, height = 4, dpi=300)
      
    ## 4 横向きグラフ - label
  g4 <- g3 + 
    geom_text(data = DataFrame, aes(x = Ratio, label = Percent_show),
              position = position_fill(vjust = 0.5), angle = 45)
      plot(g4)
          ggsave(paste0(FileName, "_4_Proportion_yoko_label.png"),
                  width= 12, height = 4, dpi=300)
  } # end function
## ------------------------------------------------------------------------------

  
  
  
## data 1 -- ざっくりした annotationで
df_ECannotxAge <- table(
  obj@meta.data$ECannot,
  obj@meta.data$Age) %>%
  data.frame() %>% 
    data.table::setnames(., c("Celltype", "Age", "CellNumber")) %>% 　#colnameを追加
    dplyr::filter(Age != "24month") %>%  
    group_by(Age) %>% 
    mutate(TotalCellNumber = sum(CellNumber)) %>% 
    ungroup() %>% 
    mutate(Ratio = round(CellNumber/TotalCellNumber, 2)) %>% 
    mutate(Percent = scales::percent_format(accuracy = 0.1)(CellNumber/TotalCellNumber)) %>% 
    # # ラベル用の割合を作る 0%はラベルからけす
    mutate(Percent_show = 
            case_when(
              .$Ratio < 0.01 ~ "",
              TRUE ~ .$Percent))
  head(df_ECannotxAge, n=20)
    write.xlsx(df_ECannotxAge, "1_CellProportion_ECannotxAge.xlsx")
## Colorの設定
Idents(obj) <- "ECannot"
(Colors_ECannot <- Idents(obj) %>% unique() %>% length() %>% hue_pal()(.) )
names(Colors_ECannot) <- levels(obj)
  Colors_ECannot

  
  
#> Execute function   
FUN.ProportionPlot(df_ECannotxAge, 
                   ClusterColors = Colors_ECannot,
                   FileName = "HeartEC_50732_ECannot",
                   GraphLabel = "ECannot")



## -------------------------------------------------------------------------------------
## ECannot_Fine ごとに表示
df_ECannot_FinexAge <- table(
  obj@meta.data$ECannot_Fine,
  obj@meta.data$Age) %>%
  data.frame() %>% 
    data.table::setnames(., c("Celltype", "Age", "CellNumber")) %>% 　#colnameを追加
    dplyr::filter(Age != "24month") %>%  
    group_by(Age) %>% 
    mutate(TotalCellNumber = sum(CellNumber)) %>% 
    ungroup() %>% 
    mutate(Ratio = round(CellNumber/TotalCellNumber, 2)) %>% 
    mutate(Percent = scales::percent_format(accuracy = 0.1)(CellNumber/TotalCellNumber)) %>% 
    # # ラベル用の割合を作る 0%はラベルからけす
    mutate(Percent_show = 
            case_when(
              .$Ratio < 0.01 ~ "",
              TRUE ~ .$Percent))
  head(df_ECannotxAge, n=20)
    write.xlsx(df_ECannotxAge, "2_CellProportion_ECannotFinexAge.xlsx")
## Colorの設定
Idents(obj) <- "ECannot_Fine"
(Colors_ECannotFine <- Idents(obj) %>% unique() %>% length() %>% hue_pal()(.) )
names(Colors_ECannotFine) <- levels(obj)
  Colors_ECannotFine
#> function   
FUN.ProportionPlot(df_ECannot_FinexAge, ClusterColors = Colors_ECannotFine,
                   FileName = "BrainEC_38540_ECannot_Fine", GraphLabel = "ECannot_Fine")


#> capillaryのみ -----------------------------------------------------------------
obj_sub <- readRDS('/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp190_AgedHeratEC_scRNAseq/Part4_EC_3,12,21month/P4-0_ECintegration_3,12,21month/HeartCap_3,12,21month_v1.rds')
  obj_sub
  # An object of class Seurat 
  # 32285 features across 35681 samples within 1 assay 
  # Active assay: RNA (32285 features, 2000 variable features)

## data frameの作成
df_Capillary_annot_FinexAge <- table(
  obj_sub@meta.data$ECannot_Fine,
  obj_sub@meta.data$Age) %>%
  data.frame() %>% 
    data.table::setnames(., c("Celltype", "Age", "CellNumber")) %>% 　#colnameを追加
    dplyr::filter(Age != "24month") %>%  
    group_by(Age) %>% 
    mutate(TotalCellNumber = sum(CellNumber)) %>% 
    ungroup() %>% 
    mutate(Ratio = round(CellNumber/TotalCellNumber, 2)) %>% 
    mutate(Percent = scales::percent_format(accuracy = 0.1)(CellNumber/TotalCellNumber)) %>% 
    # # ラベル用の割合を作る 0%はラベルからけす
    mutate(Percent_show = 
            case_when(
              .$Ratio < 0.01 ~ "",
              TRUE ~ .$Percent))
  head(df_Capillary_annot_FinexAge, n=20)
    write.xlsx(df_Capillary_annot_FinexAge, "3_Capillary_CellProportion_ECannotFinexAge.xlsx")
    
## Colorの設定
Idents(obj_sub) <- "ECannot_Fine"
(Colors_CapannotFine <- Idents(obj_sub) %>% unique() %>% length() %>% hue_pal()(.) )
names(Colors_CapannotFine) <- levels(obj_sub)
  Colors_CapannotFine
#> function   
FUN.ProportionPlot(df_Capillary_annot_FinexAge, 
                   ClusterColors = Colors_CapannotFine,
                    FileName = "HeartCapillary_35681_ECannot_Fine",
                   GraphLabel = "ECannot_Fine")