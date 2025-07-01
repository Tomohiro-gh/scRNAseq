## pySCENICのためのdataのconversion: Seruat object -> Anndata

（条件）
Seurat object (obj)　(version 5以降)を変換する

使用するのは`scCustomize`パッケージ

`scCustomize`は python環境が必要

### Conversion with scCustomize
reticulateを使って，この環境にAnndataなどconvertに必要なものをロードする


```r
library(reticulate)

print(reticulate::py_config())

```
通常、reticulateが自動で作成した環境は 'r-reticulate-env' という名前になりますが、
  -  念のため py_config() の virtualenv: のパスを確認して適切な envname を特定してください。
  -  例: /Users/tomohiro/Library/Caches/org.R-project.R/R/reticulate/uv/cache/archive-v0/cr9LxZ6Rsa0V3yD32FxUo
  -  この場合は、envname はそのパス全体を指定するか、reticulateが自動で認識する名前を使います。
  -  一度 reticulate::py_config() を実行して、virtualenv のパスを正確に確認してください。
  -  例: if virtualenv points to ".../cr9LxZ6Rsa0V3yD32FxUo", then envname should be that path or default.  
  -  重要なのは、py_config() で表示される 'python:' または 'virtualenv:' のパスが示す環境にインストールすることです。
  -  もしそれが `/Users/tomohiro/Library/Caches/org.R-project.R/R/reticulate/uv/cache/archive-v0/cr9LxZ6Rsa0V3yD32FxUo` なら、 envname にそのパスを指定するのが最も確実です。

#### この環境へ anndataをインストールする
```r
reticulate::py_install(
  packages = c("anndata", "scipy"),
  envname = "/Path/to/env/org.R-project.R/R/reticulate/uv/cache/archive-v0/cr9LxZ6Rsa0V3yD32FxUo", # py_config() で確認した仮想環境のパスを指定
  method = "auto", # reticulateに最適な方法 (conda/pip) を選ばせる
  conda = "auto",  # condaがインストールされていれば使う
  pip = TRUE,      # pipでのインストールも許可
  force = TRUE,    # 強制的に再インストールを試みる
  reinstall = TRUE # 既存のものを上書きして再インストール
)

  # 環境が正しく設定されているか確認
  print(reticulate::py_config())
```


### Check and Modify data structure
このpython環境でAnndataがインポートできているか？
```r
tryCatch({
  reticulate::py_run_string("import anndata as ad")
  message("AnnData imported successfully!")
  message(paste0("AnnData version: ", reticulate::py_eval("ad.__version__")))
}, error = function(e) {
  stop("AnnData のインポートに失敗しました。詳細: ", e$message)
})
```
AnnData imported successfully!

AnnData version: 0.11.4


```r
library(Seurat)
library(scCustomize)

wd = 'path/to/wd'

obj <- readRDS('path/to/obj.rds')
  obj

scCustomize::as.anndata(
  obj,
  file_path = wd,
  file_name = "obj.h5ad")

```
