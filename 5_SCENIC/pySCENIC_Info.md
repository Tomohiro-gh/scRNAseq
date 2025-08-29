# pySCENIC [Stein Aerts lab]

pySCENIC
#### Downloading cisTarget database
- [pySCENIC Tutorial](https://pyscenic.readthedocs.io/en/latest/installation.html)
- [GitHub](https://github.com/aertslab/pySCENIC)
- [create_cisTarget_databases](https://github.com/aertslab/create_cisTarget_databases)
- [Pypi](https://pypi.org/project/pyscenic/0.8.11/)

### Step1. cis_Target_databaseのダウンロード
マウス：　https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc_v10_clust/gene_based/

### Step2. 環境構築
[Installation](https://pyscenic.readthedocs.io/en/latest/installation.html)

Caution) 24/07/23 numpyのversionの関係で， python=3.7で作った方が良い．　python=3.10では，インストールできるものの，　`pyscenic -h`でerrorが出る

```sh
% conda create -y -n python=3.7 pyscenic
% conda activate pyscenic
% pip install pyscenic

pip install pyscenic
Collecting pyscenic
  Using cached pyscenic-0.12.1-py3-none-any.whl.metadata (9.8 kB)
Collecting ctxcore>=0.2.0 (from pyscenic)
  Using cached ctxcore-0.2.0-py3-none-any.whl.metadata (3.1 kB)
Collecting cytoolz (from pyscenic)
  Downloading cytoolz-0.12.3-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (4.6 kB)
Collecting multiprocessing-on-dill (from pyscenic)
  Using cached multiprocessing_on_dill-3.5.0a4.tar.gz (53 kB)
  Preparing metadata (setup.py) ... done
Collecting llvmlite (from pyscenic)
  Using cached llvmlite-0.39.1-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (4.7 kB)
Collecting numba>=0.51.2 (from pyscenic)
  Using cached numba-0.56.4-cp37-cp37m-macosx_10_14_x86_64.whl.metadata (2.8 kB)
Collecting attrs (from pyscenic)
  Using cached attrs-23.2.0-py3-none-any.whl.metadata (9.5 kB)
Collecting frozendict (from pyscenic)
  Using cached frozendict-2.4.4.tar.gz (315 kB)
  Installing build dependencies ... done
  Getting requirements to build wheel ... done
  Installing backend dependencies ... done
  Preparing metadata (pyproject.toml) ... done
Collecting numpy (from pyscenic)
  Using cached numpy-1.21.6-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (2.1 kB)
Collecting pandas>=1.3.5 (from pyscenic)
  Using cached pandas-1.3.5-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (12 kB)
Collecting numexpr (from pyscenic)
  Using cached numexpr-2.8.6-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (8.0 kB)
Collecting cloudpickle (from pyscenic)
  Downloading cloudpickle-2.2.1-py3-none-any.whl.metadata (6.9 kB)
Collecting dask (from pyscenic)
  Downloading dask-2022.2.0-py3-none-any.whl.metadata (3.0 kB)
Collecting distributed (from pyscenic)
  Downloading distributed-2022.2.0-py3-none-any.whl.metadata (2.9 kB)
Collecting arboreto>=0.1.6 (from pyscenic)
  Using cached arboreto-0.1.6-py2.py3-none-any.whl.metadata (5.4 kB)
Collecting boltons (from pyscenic)
  Using cached boltons-24.0.0-py3-none-any.whl.metadata (1.5 kB)
Requirement already satisfied: setuptools in ./opt/anaconda3/envs/pyscenic/lib/python3.7/site-packages (from pyscenic) (65.6.3)
Collecting pyyaml (from pyscenic)
  Using cached PyYAML-6.0.1-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (2.1 kB)
Collecting tqdm (from pyscenic)
  Using cached tqdm-4.66.4-py3-none-any.whl.metadata (57 kB)
Collecting interlap (from pyscenic)
  Using cached interlap-0.2.7.tar.gz (6.1 kB)
  Preparing metadata (setup.py) ... done
Collecting umap-learn (from pyscenic)
  Using cached umap_learn-0.5.6-py3-none-any.whl.metadata (21 kB)
Collecting loompy (from pyscenic)
  Using cached loompy-3.0.7-py3-none-any.whl
Collecting networkx (from pyscenic)
  Downloading networkx-2.6.3-py3-none-any.whl.metadata (5.0 kB)
Collecting scipy (from pyscenic)
  Using cached scipy-1.7.3-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (2.2 kB)
Collecting fsspec (from pyscenic)
  Downloading fsspec-2023.1.0-py3-none-any.whl.metadata (5.5 kB)
Collecting requests (from pyscenic)
  Using cached requests-2.31.0-py3-none-any.whl.metadata (4.6 kB)
Collecting aiohttp (from pyscenic)
  Downloading aiohttp-3.8.6-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (7.7 kB)
Collecting scikit-learn>=0.22.2 (from pyscenic)
  Using cached scikit_learn-1.0.2-cp37-cp37m-macosx_10_13_x86_64.whl.metadata (10 kB)
Collecting pyarrow>=8.0.0 (from ctxcore>=0.2.0->pyscenic)
  Downloading pyarrow-12.0.1-cp37-cp37m-macosx_10_14_x86_64.whl.metadata (3.0 kB)
Collecting importlib-metadata (from numba>=0.51.2->pyscenic)
  Using cached importlib_metadata-6.7.0-py3-none-any.whl.metadata (4.9 kB)
Collecting python-dateutil>=2.7.3 (from pandas>=1.3.5->pyscenic)
  Using cached python_dateutil-2.9.0.post0-py2.py3-none-any.whl.metadata (8.4 kB)
Collecting pytz>=2017.3 (from pandas>=1.3.5->pyscenic)
  Using cached pytz-2024.1-py2.py3-none-any.whl.metadata (22 kB)
Collecting joblib>=0.11 (from scikit-learn>=0.22.2->pyscenic)
  Using cached joblib-1.3.2-py3-none-any.whl.metadata (5.4 kB)
Collecting threadpoolctl>=2.0.0 (from scikit-learn>=0.22.2->pyscenic)
  Using cached threadpoolctl-3.1.0-py3-none-any.whl.metadata (9.2 kB)
Collecting charset-normalizer<4.0,>=2.0 (from aiohttp->pyscenic)
  Using cached charset_normalizer-3.3.2-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (33 kB)
Collecting multidict<7.0,>=4.5 (from aiohttp->pyscenic)
  Downloading multidict-6.0.5-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (4.2 kB)
Collecting async-timeout<5.0,>=4.0.0a3 (from aiohttp->pyscenic)
  Using cached async_timeout-4.0.3-py3-none-any.whl.metadata (4.2 kB)
Collecting yarl<2.0,>=1.0 (from aiohttp->pyscenic)
  Downloading yarl-1.9.4-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (31 kB)
Collecting frozenlist>=1.1.1 (from aiohttp->pyscenic)
  Downloading frozenlist-1.3.3-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (4.7 kB)
Collecting aiosignal>=1.1.2 (from aiohttp->pyscenic)
  Using cached aiosignal-1.3.1-py3-none-any.whl.metadata (4.0 kB)
Collecting asynctest==0.13.0 (from aiohttp->pyscenic)
  Downloading asynctest-0.13.0-py3-none-any.whl.metadata (5.4 kB)
Collecting typing-extensions>=3.7.4 (from aiohttp->pyscenic)
  Using cached typing_extensions-4.7.1-py3-none-any.whl.metadata (3.1 kB)
Collecting toolz>=0.8.0 (from cytoolz->pyscenic)
  Using cached toolz-0.12.1-py3-none-any.whl.metadata (5.1 kB)
Collecting packaging>=20.0 (from dask->pyscenic)
  Using cached packaging-24.0-py3-none-any.whl.metadata (3.2 kB)
Collecting partd>=0.3.10 (from dask->pyscenic)
  Using cached partd-1.4.1-py3-none-any.whl.metadata (4.6 kB)
Collecting click>=6.6 (from distributed->pyscenic)
  Using cached click-8.1.7-py3-none-any.whl.metadata (3.0 kB)
Collecting jinja2 (from distributed->pyscenic)
  Using cached jinja2-3.1.4-py3-none-any.whl.metadata (2.6 kB)
Collecting msgpack>=0.6.0 (from distributed->pyscenic)
  Downloading msgpack-1.0.5-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (8.8 kB)
Collecting psutil>=5.0 (from distributed->pyscenic)
  Using cached psutil-6.0.0-cp36-abi3-macosx_10_9_x86_64.whl.metadata (21 kB)
Collecting sortedcontainers!=2.0.0,!=2.0.1 (from distributed->pyscenic)
  Using cached sortedcontainers-2.4.0-py2.py3-none-any.whl.metadata (10 kB)
Collecting tblib>=1.6.0 (from distributed->pyscenic)
  Downloading tblib-2.0.0-py3-none-any.whl.metadata (25 kB)
Collecting zict>=0.1.3 (from distributed->pyscenic)
  Downloading zict-2.2.0-py2.py3-none-any.whl.metadata (922 bytes)
Collecting tornado>=5 (from distributed->pyscenic)
  Using cached tornado-6.2-cp37-abi3-macosx_10_9_x86_64.whl.metadata (2.5 kB)
Collecting h5py (from loompy->pyscenic)
  Using cached h5py-3.8.0-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (2.5 kB)
Collecting numpy-groupies (from loompy->pyscenic)
  Using cached numpy_groupies-0.9.22-py3-none-any.whl
Collecting dill (from multiprocessing-on-dill->pyscenic)
  Downloading dill-0.3.7-py3-none-any.whl.metadata (9.9 kB)
Collecting idna<4,>=2.5 (from requests->pyscenic)
  Using cached idna-3.7-py3-none-any.whl.metadata (9.9 kB)
Collecting urllib3<3,>=1.21.1 (from requests->pyscenic)
  Using cached urllib3-2.0.7-py3-none-any.whl.metadata (6.6 kB)
Collecting certifi>=2017.4.17 (from requests->pyscenic)
  Using cached certifi-2024.7.4-py3-none-any.whl.metadata (2.2 kB)
Collecting pynndescent>=0.5 (from umap-learn->pyscenic)
  Using cached pynndescent-0.5.13-py3-none-any.whl.metadata (6.8 kB)
Collecting locket (from partd>=0.3.10->dask->pyscenic)
  Using cached locket-1.0.0-py2.py3-none-any.whl.metadata (2.8 kB)
Collecting zipp>=0.5 (from importlib-metadata->numba>=0.51.2->pyscenic)
  Using cached zipp-3.15.0-py3-none-any.whl.metadata (3.7 kB)
Collecting six>=1.5 (from python-dateutil>=2.7.3->pandas>=1.3.5->pyscenic)
  Using cached six-1.16.0-py2.py3-none-any.whl.metadata (1.8 kB)
Collecting heapdict (from zict>=0.1.3->distributed->pyscenic)
  Downloading HeapDict-1.0.1-py3-none-any.whl.metadata (1.9 kB)
Collecting bokeh>=2.1.1 (from dask[complete]->arboreto>=0.1.6->pyscenic)
  Downloading bokeh-2.4.3-py3-none-any.whl.metadata (14 kB)
Collecting MarkupSafe>=2.0 (from jinja2->distributed->pyscenic)
  Using cached MarkupSafe-2.1.5-cp37-cp37m-macosx_10_9_x86_64.whl.metadata (3.0 kB)
Collecting pillow>=7.1.0 (from bokeh>=2.1.1->dask[complete]->arboreto>=0.1.6->pyscenic)
  Using cached Pillow-9.5.0-cp37-cp37m-macosx_10_10_x86_64.whl.metadata (9.5 kB)
Using cached pyscenic-0.12.1-py3-none-any.whl (7.1 MB)
Using cached arboreto-0.1.6-py2.py3-none-any.whl (15 kB)
Using cached ctxcore-0.2.0-py3-none-any.whl (5.8 MB)
Using cached numba-0.56.4-cp37-cp37m-macosx_10_14_x86_64.whl (2.4 MB)
Using cached llvmlite-0.39.1-cp37-cp37m-macosx_10_9_x86_64.whl (25.5 MB)
Using cached numpy-1.21.6-cp37-cp37m-macosx_10_9_x86_64.whl (16.9 MB)
Using cached pandas-1.3.5-cp37-cp37m-macosx_10_9_x86_64.whl (11.0 MB)
Using cached scikit_learn-1.0.2-cp37-cp37m-macosx_10_13_x86_64.whl (7.8 MB)
Using cached scipy-1.7.3-cp37-cp37m-macosx_10_9_x86_64.whl (33.0 MB)
Downloading aiohttp-3.8.6-cp37-cp37m-macosx_10_9_x86_64.whl (363 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 363.8/363.8 kB 21.7 MB/s eta 0:00:00
Downloading asynctest-0.13.0-py3-none-any.whl (26 kB)
Using cached attrs-23.2.0-py3-none-any.whl (60 kB)
Using cached boltons-24.0.0-py3-none-any.whl (191 kB)
Downloading cloudpickle-2.2.1-py3-none-any.whl (25 kB)
Downloading cytoolz-0.12.3-cp37-cp37m-macosx_10_9_x86_64.whl (412 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 412.1/412.1 kB 28.3 MB/s eta 0:00:00
Downloading dask-2022.2.0-py3-none-any.whl (1.1 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1.1/1.1 MB 53.2 MB/s eta 0:00:00
Downloading fsspec-2023.1.0-py3-none-any.whl (143 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 143.0/143.0 kB 13.5 MB/s eta 0:00:00
Using cached PyYAML-6.0.1-cp37-cp37m-macosx_10_9_x86_64.whl (189 kB)
Downloading distributed-2022.2.0-py3-none-any.whl (837 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 837.0/837.0 kB 41.0 MB/s eta 0:00:00
Downloading networkx-2.6.3-py3-none-any.whl (1.9 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 1.9/1.9 MB 71.2 MB/s eta 0:00:00
Using cached numexpr-2.8.6-cp37-cp37m-macosx_10_9_x86_64.whl (105 kB)
Using cached requests-2.31.0-py3-none-any.whl (62 kB)
Using cached tqdm-4.66.4-py3-none-any.whl (78 kB)
Using cached umap_learn-0.5.6-py3-none-any.whl (85 kB)
Using cached aiosignal-1.3.1-py3-none-any.whl (7.6 kB)
Using cached async_timeout-4.0.3-py3-none-any.whl (5.7 kB)
Using cached certifi-2024.7.4-py3-none-any.whl (162 kB)
Using cached charset_normalizer-3.3.2-cp37-cp37m-macosx_10_9_x86_64.whl (118 kB)
Using cached click-8.1.7-py3-none-any.whl (97 kB)
Downloading frozenlist-1.3.3-cp37-cp37m-macosx_10_9_x86_64.whl (36 kB)
Using cached idna-3.7-py3-none-any.whl (66 kB)
Using cached joblib-1.3.2-py3-none-any.whl (302 kB)
Downloading msgpack-1.0.5-cp37-cp37m-macosx_10_9_x86_64.whl (72 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 72.5/72.5 kB 5.4 MB/s eta 0:00:00
Downloading multidict-6.0.5-cp37-cp37m-macosx_10_9_x86_64.whl (28 kB)
Using cached packaging-24.0-py3-none-any.whl (53 kB)
Downloading partd-1.4.1-py3-none-any.whl (18 kB)
Using cached psutil-6.0.0-cp36-abi3-macosx_10_9_x86_64.whl (250 kB)
Downloading pyarrow-12.0.1-cp37-cp37m-macosx_10_14_x86_64.whl (24.7 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 24.7/24.7 MB 82.2 MB/s eta 0:00:00
Using cached pynndescent-0.5.13-py3-none-any.whl (56 kB)
Using cached importlib_metadata-6.7.0-py3-none-any.whl (22 kB)
Using cached python_dateutil-2.9.0.post0-py2.py3-none-any.whl (229 kB)
Using cached pytz-2024.1-py2.py3-none-any.whl (505 kB)
Using cached sortedcontainers-2.4.0-py2.py3-none-any.whl (29 kB)
Downloading tblib-2.0.0-py3-none-any.whl (11 kB)
Using cached threadpoolctl-3.1.0-py3-none-any.whl (14 kB)
Using cached toolz-0.12.1-py3-none-any.whl (56 kB)
Using cached tornado-6.2-cp37-abi3-macosx_10_9_x86_64.whl (419 kB)
Using cached typing_extensions-4.7.1-py3-none-any.whl (33 kB)
Using cached urllib3-2.0.7-py3-none-any.whl (124 kB)
Downloading yarl-1.9.4-cp37-cp37m-macosx_10_9_x86_64.whl (82 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 83.0/83.0 kB 5.9 MB/s eta 0:00:00
Downloading zict-2.2.0-py2.py3-none-any.whl (23 kB)
Downloading dill-0.3.7-py3-none-any.whl (115 kB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 115.3/115.3 kB 10.3 MB/s eta 0:00:00
Using cached h5py-3.8.0-cp37-cp37m-macosx_10_9_x86_64.whl (3.2 MB)
Using cached jinja2-3.1.4-py3-none-any.whl (133 kB)
Downloading bokeh-2.4.3-py3-none-any.whl (18.5 MB)
   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 18.5/18.5 MB 84.1 MB/s eta 0:00:00
Using cached MarkupSafe-2.1.5-cp37-cp37m-macosx_10_9_x86_64.whl (13 kB)
Using cached six-1.16.0-py2.py3-none-any.whl (11 kB)
Using cached zipp-3.15.0-py3-none-any.whl (6.8 kB)
Downloading HeapDict-1.0.1-py3-none-any.whl (3.9 kB)
Using cached locket-1.0.0-py2.py3-none-any.whl (4.4 kB)
Using cached Pillow-9.5.0-cp37-cp37m-macosx_10_10_x86_64.whl (3.4 MB)
Building wheels for collected packages: frozendict, interlap, multiprocessing-on-dill
  Building wheel for frozendict (pyproject.toml) ... done
  Created wheel for frozendict: filename=frozendict-2.4.4-cp37-cp37m-macosx_10_9_x86_64.whl size=37018 sha256=1db41ca0abba573dc99fc2ca73842957d6de33cb21c164528420a320f403aac2
  Stored in directory: /Users/tomohiro/Library/Caches/pip/wheels/51/39/85/a63fa86f9dfcd83680eb8bc9e31ce235b4bca4ca5c263d46e3
  Building wheel for interlap (setup.py) ... done
  Created wheel for interlap: filename=interlap-0.2.7-py3-none-any.whl size=6281 sha256=47b58319c31629caf4c88dfd0ce99a952e5d12159173816be89f6035649d3524
  Stored in directory: /Users/tomohiro/Library/Caches/pip/wheels/d3/6e/38/9d877cb6ae22516e2a36b01a18ded2dee5d3438b686342bca7
  Building wheel for multiprocessing-on-dill (setup.py) ... done
  Created wheel for multiprocessing-on-dill: filename=multiprocessing_on_dill-3.5.0a4-py3-none-any.whl size=64022 sha256=cf0fd7dee918ec7f4a3f66cfc5a68d2158cd6762bfa558e5bc770f788c4c31bf
  Stored in directory: /Users/tomohiro/Library/Caches/pip/wheels/ce/3c/91/f3f4e2b7aa04ca3f58d3c53cb6595f7f025c513a77e30f2bc9
Successfully built frozendict interlap multiprocessing-on-dill
Installing collected packages: sortedcontainers, pytz, msgpack, interlap, heapdict, zipp, zict, urllib3, typing-extensions, tqdm, tornado, toolz, threadpoolctl, tblib, six, pyyaml, psutil, pillow, packaging, numpy, networkx, multidict, MarkupSafe, locket, llvmlite, joblib, idna, fsspec, frozenlist, frozendict, dill, cloudpickle, charset-normalizer, certifi, boltons, asynctest, yarl, scipy, requests, python-dateutil, pyarrow, partd, numpy-groupies, numexpr, multiprocessing-on-dill, jinja2, importlib-metadata, h5py, cytoolz, async-timeout, aiosignal, scikit-learn, pandas, numba, dask, click, bokeh, attrs, pynndescent, loompy, distributed, ctxcore, aiohttp, umap-learn, arboreto, pyscenic
Successfully installed MarkupSafe-2.1.5 aiohttp-3.8.6 aiosignal-1.3.1 arboreto-0.1.6 async-timeout-4.0.3 asynctest-0.13.0 attrs-23.2.0 bokeh-2.4.3 boltons-24.0.0 certifi-2024.7.4 charset-normalizer-3.3.2 click-8.1.7 cloudpickle-2.2.1 ctxcore-0.2.0 cytoolz-0.12.3 dask-2022.2.0 dill-0.3.7 distributed-2022.2.0 frozendict-2.4.4 frozenlist-1.3.3 fsspec-2023.1.0 h5py-3.8.0 heapdict-1.0.1 idna-3.7 importlib-metadata-6.7.0 interlap-0.2.7 jinja2-3.1.4 joblib-1.3.2 llvmlite-0.39.1 locket-1.0.0 loompy-3.0.7 msgpack-1.0.5 multidict-6.0.5 multiprocessing-on-dill-3.5.0a4 networkx-2.6.3 numba-0.56.4 numexpr-2.8.6 numpy-1.21.6 numpy-groupies-0.9.22 packaging-24.0 pandas-1.3.5 partd-1.4.1 pillow-9.5.0 psutil-6.0.0 pyarrow-12.0.1 pynndescent-0.5.13 pyscenic-0.12.1 python-dateutil-2.9.0.post0 pytz-2024.1 pyyaml-6.0.1 requests-2.31.0 scikit-learn-1.0.2 scipy-1.7.3 six-1.16.0 sortedcontainers-2.4.0 tblib-2.0.0 threadpoolctl-3.1.0 toolz-0.12.1 tornado-6.2 tqdm-4.66.4 typing-extensions-4.7.1 umap-learn-0.5.6 urllib3-2.0.7 yarl-1.9.4 zict-2.2.0 zipp-3.15.0


% pyscenic -h
usage: pyscenic [-h] {grn,add_cor,ctx,aucell} ...

Single-Cell rEgulatory Network Inference and Clustering (0.12.1)

positional arguments:
  {grn,add_cor,ctx,aucell}
                        sub-command help
    grn                 Derive co-expression modules from expression matrix.
    add_cor             [Optional] Add Pearson correlations based on TF-gene
                        expression to the network adjacencies output from the
                        GRN step, and output these to a new adjacencies file.
                        This will normally be done during the "ctx" step.
    ctx                 Find enriched motifs for a gene signature and
                        optionally prune targets from this signature based on
                        cis-regulatory cues.
    aucell              Quantify activity of gene signatures across single
                        cells.

optional arguments:
  -h, --help            show this help message and exit

Arguments can be read from file using a @args.txt construct. For more
information on loom file format see http://loompy.org . For more information
on gmt file format see https://software.broadinstitute.org/cancer/software/gse
a/wiki/index.php/Data_formats .
```
