# CellBender

## Installation
Official page： https://cellbender.readthedocs.io/en/latest/installation/index.html

参考ページ：　 https://labo-code.com/bioinformatics/qc-cellbender/


--------------------
Installation on NIG supercomputer 24/4/20 -> Failed

```sh
conda create -n cellbender python=3.7
conda activate cellbender
pip install cellbender
```
しかしエラーが出てしまった．
```
ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.9' 
```

--------------------
Installation on my mac (sonoma 14.4.1) 24/04/25

```sh
conda create -n cellbender python=3.7
conda activate cellbender
pip install cellbender
```
こっちではOK. 無事runできた．

ただし，Docker imageの方は失敗した．

---------------------
##  Trial 
#### 04/25/24
下記実行
```sh
$ sh CellBenderTest_240425.sh
```

コードの中身： 

```sh
#!/bin/bash

conda activate cellbender

wd=/Volumes/../CellrangerCount/1_cellranger_count_v1
cd $wd

path_to_fastq=/Volumes/../rawfastq


ls -F $path_to_fastq | grep / | sed s/\\/$//g >> ./SampleList.txt

LOCK_FILE2="SampleList.txt"


cat ./${LOCK_FILE2} | while read line
do
    mkdir $wd/${line}_cellbender
    OutDir=$wd/${line}_cellbender

    cellbender remove-background \
        --input $wd/${line}/outs/raw_feature_bc_matrix.h5 \
        --output $OutDir/${line}_cellbender_matrix.h5
done

```
 ---------------------

### Usage
