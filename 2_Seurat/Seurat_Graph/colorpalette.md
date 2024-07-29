## Seuratで作成したcluster colorを使用する

```r
## 確認と設定
Idents(seu_obj) <- "CellTypeAnnotation"
levels(Object)

## Colorの設定
(ClusterColor <- Idents(Object) %>% unique() %>% length() %>% hue_pal()(.) )
names(ClusterColor) <- levels(Object)
  ClusterColor
```

```r
# colorの表示は 例えばclusterの数が10の場合
show(hue_pal()(10))
```

##### さらにこのcolorをbar plotなどで使用したい場合は下記のように指定すると良い
`scale_fill_manual(values = ClusterColor)` 

or

`scale_fill_manual(values = ClusterColor)` 

Example
```r
g1 <- ggplot(df, aes(x = Age, y = Ratio, fill = Celltype)) + 
    geom_bar(width = 0.6, stat = "identity", position = "fill") + 
    scale_y_continuous(expand=c(0,0), labels = percent) +
    scale_fill_manual(values = ClusterColor) +
    labs(x="", y="Ratio",title = "") + 
    theme(plot.title = element_text(hjust=0.5, size=16), # タイトルの文字サイズと位置
          plot.subtitle = element_text(size = 14, hjust = 0.5), # サブタイトルの文字サイズと位置
          axis.title.x = element_text(hjust=0.5, size= 16, face="bold"), # 軸の名前の文字の大きさ
          axis.title.y = element_text(size= 16, face="bold"),
          axis.text.x = element_text(colour="black", size= 14,  face="bold"), # 軸の目盛り文字の大きさ
          axis.text.y = element_text(colour="black", size= 14),
          axis.line=element_line(colour = "black"), # 軸の色
          axis.ticks=element_line(colour = "black"), # ティックマークの色
          panel.grid.major = element_blank(), # 主グリッドラインを描かない
          panel.grid.minor = element_blank(), # 副グリッドラインを描かない
          panel.background = element_blank(), # バックグラウンドを白にする
          text = element_text(size = 12)) +
    scale_fill_discrete(name="ECtype") 

      plot(g1)

```
