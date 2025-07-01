# pySCENIC


## setup 25/07/01　時点　----
ここが一番大事


```sh
conda create -y -n pyscenic python=3.10
conda activate pyscenic

## 特定のグレードにダウンする
pip uninstall numpy numba dask
conda install numpy=1.23.5
conda install numba=0.56.4
conda install dask=2022.2.0
``
## step3 
