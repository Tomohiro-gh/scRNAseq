# pySCENIC


## setup 25/07/01　時点　----
ポイント
#### `pip install pyscenic`を行うと，新しいパッケージが入るが，以下に不具合が出た．
- Numpy >>> `1.23.5` が良い　
- Numba >>> `0.56.4`でないとNumpyのエラーが出る
- Pandas >>> `1.5.3`がいいらしいが，`python10`でこれを使えない．現在はこのまま
- Dask >>> 現在の`2025.1`はエラーがでる
-

#### 07/01/25 ーーーー下記でうまくいった．
```sh
conda create -y -n pyscenic python=3.10
conda activate pyscenic

# まずpyscenicをPyPlでインストール
pip install pyscenic

## 特定のグレードにダウンする
pip uninstall numpy numba dask distributed
conda install numpy=1.23.5
conda install numba=0.56.4
conda install dask=2022.2.0
conda install distributed=2022.2.0

```

[ 重要！]
これでもAUCellでエラーが出る．
```
  File "/opt/anaconda3/envs/pyscenic38/lib/python3.8/site-packages/pyscenic/cli/utils.py", line 347, in append_auc_mtx

    for name, threshold in auc_thresholds.iteritems()

AttributeError: 'Series' object has no attribute 'iteritems'

```
- `iteritems`でのエラーは pandas由来．
- 非推奨的な方法だが，pyscenicのpandasの部分を書き換える >
- 参考：　[BUG]'Series' object has no attribute 'iteritems' [#475](https://github.com/aertslab/pySCENIC/issues/475#issuecomment-2489952860)

エラーのでた utils.pyをエディターで開く

  "/opt/anaconda3/envs/pyscenic_v2/lib/python3.8/site-packages/pyscenic/cli/utils.py", line 347, in append_auc_mtx

```
##書き換え前
iteritems

##
items
```

#> 最終的な環境は以下の通り．

```sh
conda list
```
#####　何度も書くが，utils.pyは書き換えている

      # packages in environment at /opt/anaconda3/envs/pyscenic_v2:
      #
      # Name                    Version                   Build  Channel
      aiohappyeyeballs          2.6.1                    pypi_0    pypi
      aiohttp                   3.12.13                  pypi_0    pypi
      aiosignal                 1.3.2                    pypi_0    pypi
      arboreto                  0.1.6                    pypi_0    pypi
      async-timeout             5.0.1                    pypi_0    pypi
      attrs                     25.3.0                   pypi_0    pypi
      bokeh                     3.7.3                    pypi_0    pypi
      boltons                   25.0.0                   pypi_0    pypi
      bzip2                     1.0.8                h99b78c6_7    conda-forge
      ca-certificates           2025.6.15            hbd8a1cb_0    conda-forge
      certifi                   2025.6.15                pypi_0    pypi
      charset-normalizer        3.4.2                    pypi_0    pypi
      click                     8.2.1                    pypi_0    pypi
      cloudpickle               3.1.1              pyhd8ed1ab_0    conda-forge
      contourpy                 1.3.2                    pypi_0    pypi
      ctxcore                   0.2.0                    pypi_0    pypi
      cytoolz                   1.0.1           py310h078409c_0    conda-forge
      dask                      2022.2.0           pyhd8ed1ab_0    conda-forge
      dask-core                 2022.2.0           pyhd8ed1ab_0    conda-forge
      dill                      0.4.0                    pypi_0    pypi
      distributed               2025.5.1                 pypi_0    pypi
      frozendict                2.4.6                    pypi_0    pypi
      frozenlist                1.7.0                    pypi_0    pypi
      fsspec                    2025.5.1           pyhd8ed1ab_0    conda-forge
      h5py                      3.14.0                   pypi_0    pypi
      idna                      3.10                     pypi_0    pypi
      importlib-metadata        8.7.0                    pypi_0    pypi
      interlap                  0.2.7                    pypi_0    pypi
      jinja2                    3.1.6              pyhd8ed1ab_0    conda-forge
      joblib                    1.5.1                    pypi_0    pypi
      lcms2                     2.17                 h7eeda09_0    conda-forge
      lerc                      4.0.0                hd64df32_1    conda-forge
      libblas                   3.9.0           32_h10e41b3_openblas    conda-forge
      libcblas                  3.9.0           32_hb3479ef_openblas    conda-forge
      libcxx                    20.1.7               ha82da77_0    conda-forge
      libdeflate                1.24                 h5773f1b_0    conda-forge
      libexpat                  2.7.0                h286801f_0    conda-forge
      libffi                    3.4.6                h1da3d7d_1    conda-forge
      libfreetype               2.13.3               hce30654_1    conda-forge
      libfreetype6              2.13.3               h1d14073_1    conda-forge
      libgfortran               5.0.0           14_2_0_h6c33f7e_103    conda-forge
      libgfortran5              14.2.0             h6c33f7e_103    conda-forge
      libjpeg-turbo             3.1.0                h5505292_0    conda-forge
      liblapack                 3.9.0           32_hc9a63f6_openblas    conda-forge
      libllvm11                 11.1.0               hfa12f05_5    conda-forge
      liblzma                   5.8.1                h39f12f2_2    conda-forge
      libopenblas               0.3.30          openmp_hf332438_0    conda-forge
      libpng                    1.6.49               h3783ad8_0    conda-forge
      libsqlite                 3.50.2               h6fb428d_0    conda-forge
      libtiff                   4.7.0                h2f21f7c_5    conda-forge
      libwebp-base              1.5.0                h2471fea_0    conda-forge
      libxcb                    1.17.0               hdb1d25a_0    conda-forge
      libzlib                   1.3.1                h8359307_2    conda-forge
      llvm-openmp               20.1.7               hdb05f8b_0    conda-forge
      llvmlite                  0.44.0                   pypi_0    pypi
      locket                    1.0.0              pyhd8ed1ab_0    conda-forge
      loompy                    3.0.8                    pypi_0    pypi
      lz4                       4.4.4                    pypi_0    pypi
      markupsafe                3.0.2           py310hc74094e_1    conda-forge
      msgpack-python            1.1.1           py310h7f4e7e6_0    conda-forge
      multidict                 6.6.3                    pypi_0    pypi
      multiprocessing-on-dill   3.5.0a4                  pypi_0    pypi
      narwhals                  1.44.0                   pypi_0    pypi
      ncurses                   6.5                  h5e97a16_3    conda-forge
      networkx                  3.4.2                    pypi_0    pypi
      numba                     0.56.4          py310h3124f1e_1    conda-forge
      numexpr                   2.11.0                   pypi_0    pypi
      numpy                     1.23.5          py310h5d7c261_0    conda-forge
      numpy-groupies            0.11.3                   pypi_0    pypi
      openjpeg                  2.5.3                h8a3d83b_0    conda-forge
      openssl                   3.5.0                h81ee809_1    conda-forge
      packaging                 25.0               pyh29332c3_1    conda-forge
      pandas                    2.3.0                    pypi_0    pypi
      partd                     1.4.2              pyhd8ed1ab_0    conda-forge
      pillow                    11.2.1          py310h61efb56_0    conda-forge
      pip                       25.1.1             pyh8b19718_0    conda-forge
      propcache                 0.3.2                    pypi_0    pypi
      psutil                    7.0.0           py310h078409c_0    conda-forge
      pthread-stubs             0.4               hd74edd7_1002    conda-forge
      pyarrow                   20.0.0                   pypi_0    pypi
      pynndescent               0.5.13                   pypi_0    pypi
      pyscenic                  0.12.1                   pypi_0    pypi
      python                    3.10.18         h6cefb37_0_cpython    conda-forge
      python-dateutil           2.9.0.post0        pyhe01879c_2    conda-forge
      python_abi                3.10                    7_cp310    conda-forge
      pytz                      2025.2             pyhd8ed1ab_0    conda-forge
      pyyaml                    6.0.2           py310hc74094e_2    conda-forge
      readline                  8.2                  h1d1bf99_2    conda-forge
      requests                  2.32.4                   pypi_0    pypi
      scikit-learn              1.7.0                    pypi_0    pypi
      scipy                     1.15.3                   pypi_0    pypi
      setuptools                59.8.0          py310hbe9552e_1    conda-forge
      six                       1.17.0             pyhd8ed1ab_0    conda-forge
      sortedcontainers          2.4.0              pyhd8ed1ab_1    conda-forge
      tblib                     3.1.0              pyhd8ed1ab_0    conda-forge
      threadpoolctl             3.6.0                    pypi_0    pypi
      tk                        8.6.13               h892fb3f_2    conda-forge
      toolz                     1.0.0              pyhd8ed1ab_1    conda-forge
      tornado                   6.5.1                    pypi_0    pypi
      tqdm                      4.67.1                   pypi_0    pypi
      typing_extensions         4.14.0             pyhe01879c_0    conda-forge
      tzdata                    2025.2                   pypi_0    pypi
      umap-learn                0.5.8                    pypi_0    pypi
      urllib3                   2.5.0                    pypi_0    pypi
      wheel                     0.45.1             pyhd8ed1ab_1    conda-forge
      xorg-libxau               1.0.12               h5505292_0    conda-forge
      xorg-libxdmcp             1.1.5                hd74edd7_0    conda-forge
      xyzservices               2025.4.0                 pypi_0    pypi
      yaml                      0.2.5                h3422bc3_2    conda-forge
      yarl                      1.20.1                   pypi_0    pypi
      zict                      3.0.0              pyhd8ed1ab_1    conda-forge
      zipp                      3.23.0                   pypi_0    pypi
      zstd                      1.5.7                h6491c7d_2    conda-forge


  
