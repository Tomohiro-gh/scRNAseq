# pySCENICのためのdataのconversion: Seruat object -> Anndata

（条件）
Seurat object (obj)　(version 5以降)を変換する

使用するのは`scCustomize`パッケージ

`scCustomize`は python環境が必要

## STEP1: Conversion with scCustomize
reticulateを使って，この環境にAnndataなどconvertに必要なものをロードする


```r
library(reticulate)

print(reticulate::py_config())

```
通常、reticulateが自動で作成した環境は 'r-reticulate-env' という名前になりますが、 念のため py_config() の virtualenv: のパスを確認して適切な envname を特定してください。
  -  例: /Users/XXXXX/Library/Caches/org.R-project.R/R/reticulate/uv/cache/archive-v0/cr9LxZ6Rsa0V3yD32FxUo
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
-------------------

## STEP2: Create `h5ad` data
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
wdへ `obj.h5ad`ファイルが保存されていれば OK


-------------------

## STEP3: Modify data structure
- 最近のトライで，このままではpyscenicが読める形ではなかった
- .loomファイルへの変換を試みる

```python
import scanpy as sc
import numpy as np
from scipy.sparse import csr_matrix, issparse
```
setting data input & output

```python

ParentDir="/Path/to/pySCENIC/"
input_h5ad_path = ParentDir + "obj.h5ad"
output_h5ad_path = ParentDir + "obj_modified.loom"

print(f"Loading AnnData from: {input_h5ad_path}")
adata = sc.read_h5ad(input_h5ad_path)
print("AnnData loaded successfully.")
```

#### adata.X が疎行列でない場合、またはCSR形式でない場合はCSRに変換 > loompyはCSR形式の疎行列を好むため
```python
if not issparse(adata.X):
    print("Converting dense adata.X to CSR matrix...")
    adata.X = csr_matrix(adata.X)
elif not isinstance(adata.X, csr_matrix):
    print("Converting sparse adata.X (non-CSR) to CSR matrix...")
    adata.X = adata.X.tocsr()
else:
    print("adata.X is already a CSR matrix. No conversion needed.")

print(f"Final adata.X type: {type(adata.X)}")
print(f"adata.X shape: {adata.X.shape}")
print(f"adata.X dtype: {adata.X.dtype}")
```

#### loompy.create() には numpy.float32 型が推奨されるため、型を変換
```python
if adata.X.dtype != np.float32:
    print("Converting adata.X dtype to float32...")
    adata.X = adata.X.astype(np.float32)
```

##### --- ここからrow_attrsとcol_attrsの定義 ---
- AnnDataのobs (細胞メタデータ) と var (遺伝子メタデータ) をDataFrameに変換
- Loomファイルには `obs` が `col_attrs`、`var` が `row_attrs` に対応します
- 必ずDataFrame形式である必要があります


```python
# adata.obs (細胞属性) -> col_attrs (Loomファイルの列属性)
col_attrs = {}
for col in adata.obs.columns:
    # CategoricalDtypeの場合、文字列に変換してNumPy配列にする
    if pd.api.types.is_categorical_dtype(adata.obs[col]):
        print(f"Converting categorical column '{col}' in adata.obs to string array.")
        col_attrs[col] = np.array(adata.obs[col].astype(str).values)
    else:
        # それ以外のデータ型はこれまで通り
        col_attrs[col] = np.array(adata.obs[col].values, 
                                   dtype=str if adata.obs[col].dtype == 'object' else adata.obs[col].dtype)
col_attrs['CellID'] = np.array(adata.obs_names.values, dtype=str)

```

```python
# adata.var (遺伝子属性) -> row_attrs (Loomファイルの行属性)
row_attrs = {}
for col in adata.var.columns:
    # CategoricalDtypeの場合、文字列に変換してNumPy配列にする
    if pd.api.types.is_categorical_dtype(adata.var[col]):
        print(f"Converting categorical column '{col}' in adata.var to string array.")
        row_attrs[col] = np.array(adata.var[col].astype(str).values)
    else:
        # それ以外のデータ型はこれまで通り
        row_attrs[col] = np.array(adata.var[col].values, 
                                   dtype=str if adata.var[col].dtype == 'object' else adata.var[col].dtype)
row_attrs['Gene'] = np.array(adata.var_names.values, dtype=str)
```

##### --- Loomファイルへの書き出し ---
```python
print(f"Writing to Loom file: {output_loom_path}")
try:
    # loompy.create の行列は (genes x cells) 形式を期待するため、adata.X を転置 (adata.X は (cells x genes))
    # attrsは辞書のリストではなく、辞書のNumPy arrayを値として持つ形が最も堅牢です
    loompy.create(output_loom_path, adata.X.T, row_attrs=row_attrs, col_attrs=col_attrs)
    print("Loom file created successfully.")
except Exception as e:
    print(f"Error creating Loom file: {e}")
    print("Please review the structure of adata.obs and adata.var for compatibility with loompy.")
    print("Ensure all attributes are convertible to NumPy arrays.")

```



    

