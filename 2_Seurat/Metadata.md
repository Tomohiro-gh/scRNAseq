## Handling Metadata

#### Metadataの追加
- `AddMetaData`関数
```r
# Example1
SampleName <- "AdultZebrafish"
obj <- AddMetaData(object = obj, metadata = SampleName, col.name = "SampleName")

# Example2

```


#### 既存のcolumnの情報を参考に新しいカラムへ書き換える
- 　dplyr packageの `mutate`と`case_when`を使用する
- 　`TRUE`で上記で指定したもの以外全てをさす
```r
# Example1
ECs@meta.data <- ECs@meta.data %>% 
  mutate(dpf = case_when(.$hpf %in% c(16, 18, 21) ~ "0.66 - 0.87",
                   .$hpf %in% c(24, 26, 28, 30) ~ "1.0 - 1.25",
                   .$hpf %in% c(32, 34, 36, 38, 40) ~ "1.33 - 1.66",
                   .$hpf %in% c(42, 44, 46, 48) ~ "1.75 - 2.0",
                   .$hpf %in% c(50, 52, 54, 56, 58, 60) ~ "2.08 - 2.5",
                   .$hpf %in% c(62, 64, 66, 68, 70, 72) ~ "2.58 - 3.0",
                   .$hpf %in% c(74, 76, 78, 80, 82, 84) ~ "3.08 - 3.5",
                   .$hpf %in% c(86, 88, 90, 92, 94, 96) ~ "3.58 - 4.0",
                   TRUE ~ "4.0 - 5.0"))
```
