# Velocyto


## Installation

Install error at 24/04/20 
```
Traceback (most recent call last):
  File "/home/tomohiro/anaconda3/envs/velo/bin/velocyto", line 5, in <module>
    from velocyto.commands.velocyto import cli
  File "/home/tomohiro/anaconda3/envs/velo/lib/python3.9/site-packages/velocyto/__init__.py", line 13, in <module>
    from .estimation import fit_slope, _fit1_slope, clusters_stats
  File "/home/tomohiro/anaconda3/envs/velo/lib/python3.9/site-packages/velocyto/estimation.py", line 7, in <module>
    from .speedboosted import _colDeltaCor, _colDeltaCorLog10, _colDeltaCorSqrt
ImportError: /home/tomohiro/anaconda3/envs/velo/lib/python3.9/site-packages/velocyto/speedboosted.cpython-39-x86_64-linux-gnu.so: undefined symbol: __log10_finite
```
