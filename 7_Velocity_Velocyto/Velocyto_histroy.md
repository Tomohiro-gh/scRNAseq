# Velocytoの環境構築記録

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


