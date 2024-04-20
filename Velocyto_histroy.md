## Velocytoに関する環境構築の記録

velocytoのための環境構築の備忘録

----
### NIG supercomputer on 04/20/24
なぜかうごかくなくなった． velocytoしたらこれ．
```
$ velocyto


```
そこで新たに環境構築し直してみる



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
