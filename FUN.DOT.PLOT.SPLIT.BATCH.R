
message(
"FUN.DOT.PLOT.SPLIT.BATCH(
  obj,
  features, 
  Group_oi, 
  Split_oi, 
  grad_color = NULL,
  fname_body,
  w, 
  h
)")
#> 3. 遺伝子の発現量と発現割合を示す dotplot (Seurat package)
#> seuratのdotplotパッケージを使って，データを取得し，ggplotで自分で描写する
FUN.DOT.PLOT.SPLIT.BATCH <- function(obj,
                                     features, 
                                     Group_oi, 
                                     Split_oi, 
                                     grad_color = NULL,
                                     fname_body,
                                     w, h){
  require(Seurat)
  require(tidyverse)
  
  #> objectで Group_oi x Split_oiのカラムを作る
  vals <- obj |>
    FetchData(vars = c(Group_oi, Split_oi)) #|>
    # set_names(c("GROUP_oi", "SPLIT_oi"))
  #> Group_oixSplit_oiのカラムを作る
  obj$tmp <- paste(vals[[Group_oi]], vals[[Split_oi]], sep = "_x_")
  print(head(obj$tmp))
  
  ##> ここからloop --
  for (a_gene in features){ 
    # seuratのdot plotでdataを取得
    dp <- 
      obj |> 
      DotPlot(a_gene,
              group.by = "tmp",
              scale = FALSE)

    # 平均発現量とパーセントを取り出す
    df_tmp <-
      dp@data |> 
      dplyr::mutate(
        GROUP_oi = str_split(pattern = "_x_", id, simplify = TRUE)[,1],
        SPLIT_oi = str_split(pattern = "_x_", id, simplify = TRUE)[,2]) |>
      dplyr::mutate(
        GROUP_oi = forcats::fct_relevel(
          factor(GROUP_oi), 
          levels(obj@meta.data[[Group_oi]])), #dotで表示するときは逆にする
        SPLIT_oi = forcats::fct_relevel(
          factor(SPLIT_oi), 
          rev(levels(obj@meta.data[[Split_oi]]))))

    #> csvで保存
    # df_tmp |>
    #   # dplyr::select(avg.exp, avg.exp.scaled, pct.exp, Type, Age) |>
    #   dplyr::select(avg.exp, pct.exp, Type, Age) |>
    #   # Ageを列に展開>  avg.exp と pct.exp の両方を使う
    #   pivot_wider(
    #     names_from = Age, 
    #     values_from = c(avg.exp, pct.exp),
    #     names_glue = "{Age}_{.value}"
    #   )|>
    #   #> 並び替え
    #   dplyr::select(
    #     Type,
    #     # Age順に並び替え
    #     dplyr::all_of(as.vector(
    #       t(outer(levels(obj$Age), c("_avg.exp", "_pct.exp"),
    #               paste0)))) ) |>
    #   # 列名を「Young_avg.exp」のような形式にする
    #   mutate(Type = forcats::fct_relevel(
    #     factor(Type), levels(obj$Type))) |>
    #   arrange(Type) |>
    #   #小数点4けたまでに
    #   dplyr::mutate(across(where(is.numeric), \(x) round(x, 4))) |>
    #   #> save
    #   write.csv(paste0("3_dotplot_", fname_body,"_", a_gene,".csv"),
    #             row.names = FALSE)
    
    #> graph - ggplot　
    #> dot plot color
    if(is.null(grad_color)==TRUE){
      grad_color = colorRampPalette(colors = c("grey", "blue"))(20)
    }
    #> 
    df_tmp |>
      ggplot(aes(x = GROUP_oi, y = SPLIT_oi)) + 
      geom_point(aes(size = pct.exp, color = avg.exp)) +
      scale_color_gradientn(colours = grad_color) + # gradient color
      scale_size(range = c(2, 8)) +
      scale_x_discrete(position = "top") + 
      labs(title = a_gene,
           x = "",  y = "", 
           size = "pct.exp.",
           color = "ave.exp.") +
      guides(
        size = guide_legend(order = 1),  # 1番目に表示（上）
        color = guide_colorbar(order = 2) # 2番目に表示（下）
      ) +
      theme_minimal() +
      theme(
        title = element_text(size = 18, color = "black", hjust = 1, vjust = 1),
        axis.text = element_text(size = 18,color = "black"),
        axis.text.x = element_text(
          hjust = 0, vjust = 0, angle = 45),
        axis.text.y = element_text(
          hjust = 1),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        # .legend ...
        legend.key.size = unit(0.4, "cm"),   # 凡例のアイコン全体のサイズ
        legend.key.height = unit(0.4, "cm"), # 高さを低くしたい場合
        legend.key.width = unit(0.4, "cm"),   # 凡例の幅を狭くしたい場合
        # legend.position = "top",
        panel.grid = element_blank(),
        panel.background = element_blank(),
        aspect.ratio = 0.5) 
    ggsave(paste0("DotPlot_", fname_body,"_",a_gene, ".png"),
           width = w, height = h, dpi = 300)
  } # end loop 
} ## end function ---

#> 実行
# FUN.DOT.PLOT.BATCH(t5, test_genes, "alltissues")
# FUN.DOT.PLOT.BATCH(ht, c("Aplnr"), "Heart")
# FUN.DOT.PLOT.BATCH(br, c("Smad7"), "Brain")
# FUN.DOT.PLOT.BATCH(bec, c("Smad7"), "Brain-Capillary")

message("関数をロードしました")
message("example:")
message("my_colors <- 
        scale_color_gradientn(
          colours = wes_palette('Zissou1', type = 'continuous'))")
message("HeartCap |> 
          FUN.DOT.PLOT.SPLIT.BATCH(
            features,
            'ECannot_Fine',
            'Age',
            grad_color = my_colors,
            'HeartCaps'
            8, 5)")


