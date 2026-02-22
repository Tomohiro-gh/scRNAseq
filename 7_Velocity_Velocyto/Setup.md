# Velocytoの環境構築記録


----
### 02/22/26 on macbook (M3)
- ポイントは`python=3.8`
以下で解決．

#### 26-02-22 以下でOK
```sh
conda create -n velocyto python=3.8
conda activate velocyto

# preparation before installing velocyto
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install pysam

#> この後に，
pip install velocyto
```

確認
```sh
% velocyto --help     
OMP: Info #276: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead.
Usage: velocyto [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  run            Runs the velocity analysis outputting a loom file
  run10x         Runs the velocity analysis for a Chromium Sample
  run-dropest    Runs the velocity analysis on DropEst preprocessed data
  run-smartseq2  Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)
  tools          helper tools for velocyto


% conda list                                                               
# packages in environment at /opt/anaconda3/envs/celloracle:
#
# Name                    Version                   Build  Channel
_openmp_mutex             4.5                  7_kmp_llvm    conda-forge
brotli                    1.2.0                h7d5ae5b_1    conda-forge
brotli-bin                1.2.0                hc919400_1    conda-forge
brotli-python             1.0.9            py38h2b1e499_8    conda-forge
bzip2                     1.0.8                hd037594_9    conda-forge
c-ares                    1.34.6               hc919400_0    conda-forge
ca-certificates           2026.1.4             hbd8a1cb_0    conda-forge
cached-property           1.5.2                hd8ed1ab_1    conda-forge
cached_property           1.5.2              pyha770c72_1    conda-forge
certifi                   2024.8.30          pyhd8ed1ab_0    conda-forge
cffi                      1.17.0           py38h858044d_0    conda-forge
charset-normalizer        3.4.0              pyhd8ed1ab_0    conda-forge
click                     8.1.7           unix_pyh707e725_0    conda-forge
contourpy                 1.1.1            py38h9afee92_1    conda-forge
cycler                    0.12.1             pyhd8ed1ab_0    conda-forge
cython                    3.0.11           py38h11842c7_0    conda-forge
fonttools                 4.53.1           py38h3237794_0    conda-forge
freetype                  2.14.1               hce30654_0    conda-forge
h2                        4.1.0              pyhd8ed1ab_0    conda-forge
h5py                      3.11.0          nompi_py38h70800c8_102    conda-forge
hdf5                      1.14.3          nompi_ha698983_109    conda-forge
hpack                     4.0.0              pyh9f0ad1d_0    conda-forge
hyperframe                6.0.1              pyhd8ed1ab_0    conda-forge
icu                       78.2                 h38cb7af_0    conda-forge
idna                      3.10               pyhd8ed1ab_0    conda-forge
importlib-metadata        8.5.0              pyha770c72_0    conda-forge
importlib-resources       6.4.5              pyhd8ed1ab_0    conda-forge
importlib_resources       6.4.5              pyhd8ed1ab_0    conda-forge
joblib                    1.4.2              pyhd8ed1ab_0    conda-forge
kiwisolver                1.4.5            py38h9afee92_1    conda-forge
krb5                      1.22.2               h385eeb1_0    conda-forge
lcms2                     2.18                 hdfa7624_0    conda-forge
lerc                      4.0.0                hd64df32_1    conda-forge
libaec                    1.1.5                h8664d51_0    conda-forge
libblas                   3.9.0           20_osxarm64_openblas    conda-forge
libbrotlicommon           1.2.0                hc919400_1    conda-forge
libbrotlidec              1.2.0                hc919400_1    conda-forge
libbrotlienc              1.2.0                hc919400_1    conda-forge
libcblas                  3.9.0           20_osxarm64_openblas    conda-forge
libcurl                   8.18.0               hd5a2499_1    conda-forge
libcxx                    21.1.8               h55c6f16_2    conda-forge
libdeflate                1.25                 hc11a715_0    conda-forge
libedit                   3.1.20250104    pl5321hafb1f1b_0    conda-forge
libev                     4.33                 h93a5062_2    conda-forge
libffi                    3.5.2                hcf2aa1b_0    conda-forge
libfreetype               2.14.1               hce30654_0    conda-forge
libfreetype6              2.14.1               h6da58f4_0    conda-forge
libgcc                    15.2.0              hcbb3090_18    conda-forge
libgfortran               15.2.0              h07b0088_18    conda-forge
libgfortran5              15.2.0              hdae7583_18    conda-forge
libjpeg-turbo             3.1.2                hc919400_0    conda-forge
liblapack                 3.9.0           20_osxarm64_openblas    conda-forge
libllvm14                 14.0.6               hd1a9a77_4    conda-forge
liblzma                   5.8.2                h8088a28_0    conda-forge
liblzma-devel             5.8.2                h8088a28_0    conda-forge
libnghttp2                1.67.0               hc438710_0    conda-forge
libopenblas               0.3.25          openmp_h6c19121_0    conda-forge
libpng                    1.6.55               h132b30e_0    conda-forge
libsqlite                 3.51.2               h1ae2325_0    conda-forge
libssh2                   1.11.1               h1590b86_0    conda-forge
libtiff                   4.7.1                h4030677_1    conda-forge
libwebp-base              1.6.0                h07db88b_0    conda-forge
libxcb                    1.17.0               hdb1d25a_0    conda-forge
libzlib                   1.3.1                h8359307_2    conda-forge
llvm-openmp               21.1.8               h4a912ad_0    conda-forge
llvmlite                  0.41.1           py38hbed3d3f_0    conda-forge
loompy                    3.0.8                    pypi_0    pypi
matplotlib                3.7.3            py38h150bfb4_0    conda-forge
matplotlib-base           3.7.3            py38hef9d0d7_0    conda-forge
munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
ncurses                   6.5                  h5e97a16_3    conda-forge
numba                     0.58.1           py38hbe32c8a_0    conda-forge
numpy                     1.24.4           py38ha84db1f_0    conda-forge
numpy-groupies            0.9.22                   pypi_0    pypi
openjpeg                  2.5.4                hbfb3c88_0    conda-forge
openssl                   3.6.1                hd24854e_1    conda-forge
packaging                 26.0               pyhcf101f3_0    conda-forge
pandas                    2.0.3                    pypi_0    pypi
pillow                    10.4.0           py38h2c6aaed_0    conda-forge
pip                       24.3.1             pyh8b19718_0    conda-forge
platformdirs              4.3.6              pyhd8ed1ab_0    conda-forge
pooch                     1.8.2              pyhd8ed1ab_0    conda-forge
pthread-stubs             0.4               hd74edd7_1002    conda-forge
pycparser                 2.22               pyhd8ed1ab_0    conda-forge
pyparsing                 3.1.4              pyhd8ed1ab_0    conda-forge
pysam                     0.23.3                   pypi_0    pypi
pysocks                   1.7.1              pyha2e5f31_6    conda-forge
python                    3.8.20          h7d35d02_2_cpython    conda-forge
python-dateutil           2.9.0              pyhd8ed1ab_0    conda-forge
python_abi                3.8                      8_cp38    conda-forge
pytz                      2025.2                   pypi_0    pypi
readline                  8.3                  h46df422_0    conda-forge
requests                  2.32.3             pyhd8ed1ab_0    conda-forge
scikit-learn              1.3.2            py38he1bc1c9_2    conda-forge
scipy                     1.10.1           py38h038e806_3    conda-forge
setuptools                75.3.0             pyhd8ed1ab_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
threadpoolctl             3.5.0              pyhc1e730c_0    conda-forge
tk                        8.6.13               h010d191_3    conda-forge
tornado                   6.4.1            py38h3237794_0    conda-forge
tzdata                    2025.3                   pypi_0    pypi
unicodedata2              15.1.0           py38hb192615_0    conda-forge
urllib3                   2.2.3              pyhd8ed1ab_0    conda-forge
velocyto                  0.17.17                  pypi_0    pypi
wheel                     0.45.1             pyhd8ed1ab_0    conda-forge
xorg-libxau               1.0.12               hc919400_1    conda-forge
xorg-libxdmcp             1.1.5                hc919400_1    conda-forge
xz                        5.8.2                hd0f0c4f_0    conda-forge
xz-gpl-tools              5.8.2                hd0f0c4f_0    conda-forge
xz-tools                  5.8.2                h8088a28_0    conda-forge
zipp                      3.21.0             pyhd8ed1ab_0    conda-forge
zstandard                 0.19.0           py38hb991d35_0    conda-forge
zstd                      1.5.7                hbf9d68e_6    conda-forge

```


----
### 08/16/24 on macstudio <- inprogress

```sh
conda create -n volo3.6 python=3.6
conda create -n volo3.6 python=3.6


conda install numpy
numpy scipy cython numba matplotlib scikit-learn h5py click
conda 
```

----
### 04/20/24 on NIG supercomputer 
なぜかうごかくなくなった． velocytoしたらこれ． -> 24/04/07では使えていたのになぜ？？？
```sh
$ velocyto
## error message
Traceback (most recent call last):
  File "/home/tomohiro/anaconda3/envs/velo3.8/bin/velocyto", line 7, in <module>
    from velocyto.commands.velocyto import cli
  File "/home/tomohiro/anaconda3/envs/velo3.8/lib/python3.8/site-packages/velocyto/__init__.py", line 12, in <module>
    from .neighbors import BalancedKNN, convolve_by_sparse_weights
  File "/home/tomohiro/anaconda3/envs/velo3.8/lib/python3.8/site-packages/velocyto/neighbors.py", line 3, in <module>
    from sklearn.neighbors import kneighbors_graph, NearestNeighbors
  File "/home/tomohiro/.local/lib/python3.8/site-packages/sklearn/__init__.py", line 80, in <module>
    from .base import clone
  File "/home/tomohiro/.local/lib/python3.8/site-packages/sklearn/base.py", line 21, in <module>
    from .utils import _IS_32BIT
  File "/home/tomohiro/.local/lib/python3.8/site-packages/sklearn/utils/__init__.py", line 20, in <module>
    from scipy.sparse import issparse
  File "/home/tomohiro/.local/lib/python3.8/site-packages/scipy/sparse/__init__.py", line 229, in <module>
    from .base import *
  File "/home/tomohiro/.local/lib/python3.8/site-packages/scipy/sparse/base.py", line 8, in <module>
    from .sputils import (isdense, isscalarlike, isintlike,
  File "/home/tomohiro/.local/lib/python3.8/site-packages/scipy/sparse/sputils.py", line 17, in <module>
    supported_dtypes = [np.typeDict[x] for x in supported_dtypes]
  File "/home/tomohiro/.local/lib/python3.8/site-packages/scipy/sparse/sputils.py", line 17, in <listcomp>
    supported_dtypes = [np.typeDict[x] for x in supported_dtypes]
  File "/home/tomohiro/anaconda3/envs/velo3.8/lib/python3.8/site-packages/numpy/__init__.py", line 320, in __getattr__
    raise AttributeError("module {!r} has no attribute "
AttributeError: module 'numpy' has no attribute 'typeDict'

```
そこで新たに環境構築し直してみる

```sh
conda create -n velo3.6 python=3.6
conda activate velo3.6

conda install numpy scipy cython numba matplotlib scikit-learn h5py click
pip install pysam

pip install -U --no-deps velocyto
```

確認 -> エラー

```sh
$ velocyto --help
## (error)
Traceback (most recent call last):                                                                                                                                                
  File "/home/tomohiro/anaconda3/envs/velo3.6/bin/velocyto", line 5, in <module>                                                                                                  
    from velocyto.commands.velocyto import cli                                                                                                                                    
  File "/home/tomohiro/anaconda3/envs/velo3.6/lib/python3.6/site-packages/velocyto/__init__.py", line 15, in <module>                                                             
    from .analysis import VelocytoLoom, scatter_viz, ixs_thatsort_a2b, load_velocyto_hdf5                                                                                         
  File "/home/tomohiro/anaconda3/envs/velo3.6/lib/python3.6/site-packages/velocyto/analysis.py", line 16, in <module>                                                             
    import loompy                                                                                                                                                                 
ModuleNotFoundError: No module named 'loompy'
###                                                                             
```
loompyを入れたところ問題なかった

```sh
$ conda install loompy
$ velocyto --help

# OK!!!! no error
# Usage: velocyto [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  run            Runs the velocity analysis outputting a loom file
  run10x         Runs the velocity analysis for a Chromium Sample
  run-dropest    Runs the velocity analysis on DropEst preprocessed data
  run-smartseq2  Runs the velocity analysis on SmartSeq2 data (independent bam file per cell)
  tools          helper tools for velocyto
```
-> OK! 無事にrunが確認できた (Exp190-1and-2_HeartEC_230419.sh)

----

----
### NIG supercomputer on 01/19/24

python version 3.8.18

velocyto :biocondaでインストール
conda install velocyto

```sh
conda create -n volo3.8.18 python=3.8.18
conda activate volo3.8.18
## 入れる前に
pip install pysam
conda install numpy scipy cython numba matplotlib scikit-learn h5py click



```

----


----
### NIG supercomputer on 01/19/24
velocyto installation guide -> http://velocyto.org/velocyto.py/install/index.html
```sh
conda create -n velo python=3.9.16
conda activate velo

pip install pysam
conda install numpy scipy cython numba matplotlib scikit-learn h5py click
```
##### install velocyto
```sh
pip install velocyto
# or
pip install -U --no-deps velocyto  
```
Tips: velocyto 0.17 is an alpha release and it is updated often. If you installed with pip make sure you run pip install -U --no-deps velocyto now and then.

Install with conda -> This installation method is NOT currently available. Our plan is make it available upon the 1.0 release.

インストールエラーに関する記述 from velocyto.commands.velocyto import cli

velocyto-team/velocyto.py#53

velocyto-team/velocyto.py#186

----

--------

### 環境構築後のvelocytoの準備 ~ 2022 year

#### Step1. mappingに使用したgtfファイルを用意する
gtfファイルをworking directoryへコピーしておく

#### Step2. masked gtfの準備
繰り返し配列をはじめとした特定の領域をマスクしたgtf ,  genome dataを得る
-> velocytoではこの情報を使用することを推奨している

これを作るためには[repeat masker](http://www.repeatmasker.org/)というものが有名らしい -> [取得の仕方](https://labs.wsu.edu/winuthayanon/scrna/how-to-analyze-single%e2%80%90cell-rna%e2%80%90seq/explaining-velocyto-command-line/
)

##### 実際の実行 (zebrafishを例に)
1. [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=611454127_NtvlaW6xBSIRYJEBI0iRDEWisITa&clade=mammal&org=Mouse&db=mm10&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56694976-56714605&hgta_outputType=primaryTable&hgta_outputType=gff&hgta_outFileName=mm10_rmsk.gtf) へアクセス　 

2. optionの設定
   - clade : Vertebrate
   - genome : Zebrafish を指定
   - assmebly: GRCz11/DanRer11
   - Track: RepeatMasker
   - tableで　rmsk
   - region : genome
   - output format : GTF
   - output filename : 名前を指定

3. get outputで作成


