# CellBender

## Installation
Official page： https://cellbender.readthedocs.io/en/latest/installation/index.html

Japanese：　 https://labo-code.com/bioinformatics/qc-cellbender/

10X: https://www.10xgenomics.com/jp/analysis-guides/background-removal-guidance-for-single-cell-gene-expression-datasets-using-third-party-tools


### 参考ページ
https://rpubs.com/kenneditodd/doublet_finder_example


--------------------
#### Installation on NIG supercomputer
>  24/4/20 -> Failed
> しかしエラーが出てしまった．
> ```
> ImportError: /lib64/libstdc++.so.6: version `CXXABI_1.3.9' 
> ```
> pending

Update 24/04/26 : mac studioと同じようにしたら無事install完了

```sh
conda create -n CellBender python=3.7
conda activate CellBender
pip install cellbender
```


-------------------- 
#### Installation on my mac (sonoma 14.4.1)
24/04/25 -> suceeded!

```sh
conda create -n CellBender python=3.7
conda activate CellBender
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
### Tutorial
#### Example1
 Expected cellのパラメータのみ設定
```sh
conda activate cellbender

SampleName=MySampleName

cellbender remove-background \
        --input $InputDir/raw_feature_bc_matrix.h5 \
        --output $OutDir/$SampleName_cellbender_matrix.h5 \
        --expected-cells 9500
```

#### Example2
・　Dropletの数を指定
```sh
conda activate cellbender

SampleName=MySampleName

cellbender remove-background \
        --input $InputDir/raw_feature_bc_matrix.h5 \
        --output $OutDir/$SampleName_cellbender_matrix.h5 \
        --expected-cells 9500 \
        --total-droplets-included 30000 \
```    

---------------------

### 出力ファイル
・　${SampleName}_cellbender_matrix.h5 : 正規化された環境トランスクリプトの豊富さ、各ドロップレットの汚染度合い、バックグラウンド補正後の遺伝子発現の低次元埋め込み、およびバックグラウンド補正後のカウント行列（CSC疎行列形式）などが含まれるraw dataのファイル

・　${SampleName}_cellbender_matrix_filtered.h5 : 後方の細胞確率が0.5を超えるドロップレットを含む -> Seuratで解析するファイルだが，現在のcelranger version 7のフォーマットでは受け付けてくれない．

実は現在のsueratのversionでは`Read10X_h5()`関数で読めないらしい -> [Documentation](https://cellbender.readthedocs.io/en/latest/tutorial/index.html#open-in-seurat)

> Seurat 4.0.2 uses a dataloader Read10X_h5() which is not currently compatible with the CellBender output file format.

[Issues#315](https://github.com/broadinstitute/CellBender/issues/315) も参照

解決方法がTutorialに書いてある： 

```sh
$ ptrepack --complevel 5 tiny_output_filtered.h5:/matrix tiny_output_filtered_seurat.h5:/matrix
```
これを通せば，`Read10X_h5()`で読めるようになる．


---------------------
### Usage

・　Input file: raw_feature_bc_matrix.h5　（cellranger countの出力ファイル) <- これのみでいい


```

$ cellbender remove-background -h

usage: cellbender remove-background [-h] --input INPUT_FILE --output
                                    OUTPUT_FILE [--cuda]
                                    [--checkpoint INPUT_CHECKPOINT_TARBALL]
                                    [--expected-cells EXPECTED_CELL_COUNT]
                                    [--total-droplets-included TOTAL_DROPLETS]
                                    [--force-cell-umi-prior FORCE_CELL_UMI_PRIOR]
                                    [--force-empty-umi-prior FORCE_EMPTY_UMI_PRIOR]
                                    [--model {naive,simple,ambient,swapping,full}]
                                    [--epochs EPOCHS]
                                    [--low-count-threshold LOW_COUNT_THRESHOLD]
                                    [--z-dim Z_DIM]
                                    [--z-layers Z_HIDDEN_DIMS [Z_HIDDEN_DIMS ...]]
                                    [--training-fraction TRAINING_FRACTION]
                                    [--empty-drop-training-fraction FRACTION_EMPTIES]
                                    [--ignore-features BLACKLISTED_GENES [BLACKLISTED_GENES ...]]
                                    [--fpr FPR [FPR ...]]
                                    [--exclude-feature-types EXCLUDE_FEATURES [EXCLUDE_FEATURES ...]]
                                    [--projected-ambient-count-threshold AMBIENT_COUNTS_IN_CELLS_LOW_LIMIT]
                                    [--learning-rate LEARNING_RATE]
                                    [--checkpoint-mins CHECKPOINT_MIN]
                                    [--final-elbo-fail-fraction FINAL_ELBO_FAIL_FRACTION]
                                    [--epoch-elbo-fail-fraction EPOCH_ELBO_FAIL_FRACTION]
                                    [--num-training-tries NUM_TRAINING_TRIES]
                                    [--learning-rate-retry-mult LEARNING_RATE_RETRY_MULT]
                                    [--posterior-batch-size POSTERIOR_BATCH_SIZE]
                                    [--posterior-regularization {PRq,PRmu,PRmu_gene}]
                                    [--alpha PRQ_ALPHA] [--q CDF_THRESHOLD_Q]
                                    [--estimator {map,mean,cdf,sample,mckp}]
                                    [--estimator-multiple-cpu]
                                    [--constant-learning-rate]
                                    [--cpu-threads N_THREADS] [--debug]
                                    [--truth TRUTH_FILE]

Remove background RNA from count matrix.

optional arguments:
  -h, --help            show this help message and exit
  --input INPUT_FILE    Data file on which to run tool. Data must be un-
                        filtered: it should include empty droplets. The
                        following input formats are supported: CellRanger v2
                        and v3 (.h5 or the directory that contains the .mtx
                        file), Dropseq DGE (.txt or .txt.gz), BD Rhapsody
                        (.csv or .csv.gz), AnnData (.h5ad), and Loom (.loom).
                        (default: None)
  --output OUTPUT_FILE  Output file location (the path must exist, and the
                        file name must have .h5 extension). (default: None)
  --cuda                Including the flag --cuda will run the inference on a
                        GPU. (default: False)
  --checkpoint INPUT_CHECKPOINT_TARBALL
                        Checkpoint tarball produced by the same version of
                        CellBender remove-background. If present, and the
                        workflow hashes match, training will restart from this
                        checkpoint. (default: ckpt.tar.gz)
  --expected-cells EXPECTED_CELL_COUNT
                        Number of cells expected in the dataset (a rough
                        estimate within a factor of 2 is sufficient).
                        (default: None)
  --total-droplets-included TOTAL_DROPLETS
                        The number of droplets from the rank-ordered UMI plot
                        that will have their cell probabilities inferred as an
                        output. Include the droplets which might contain
                        cells. Droplets beyond TOTAL_DROPLETS_INCLUDED should
                        be 'surely empty' droplets. (default: None)
  --force-cell-umi-prior FORCE_CELL_UMI_PRIOR
                        Ignore CellBender's heuristic prior estimation, and
                        use this prior for UMI counts in cells. (default:
                        None)
  --force-empty-umi-prior FORCE_EMPTY_UMI_PRIOR
                        Ignore CellBender's heuristic prior estimation, and
                        use this prior for UMI counts in empty droplets.
                        (default: None)
  --model {naive,simple,ambient,swapping,full}
                        Which model is being used for count data. 'naive'
                        subtracts the estimated ambient profile. 'simple' does
                        not model either ambient RNA or random barcode
                        swapping (for debugging purposes -- not recommended).
                        'ambient' assumes background RNA is incorporated into
                        droplets. 'swapping' assumes background RNA comes from
                        random barcode swapping (via PCR chimeras). 'full'
                        uses a combined ambient and swapping model. (default:
                        full)
  --epochs EPOCHS       Number of epochs to train. (default: 150)
  --low-count-threshold LOW_COUNT_THRESHOLD
                        Droplets with UMI counts below this number are
                        completely excluded from the analysis. This can help
                        identify the correct prior for empty droplet counts in
                        the rare case where empty counts are extremely high
                        (over 200). (default: 5)
  --z-dim Z_DIM         Dimension of latent variable z. (default: 64)
  --z-layers Z_HIDDEN_DIMS [Z_HIDDEN_DIMS ...]
                        Dimension of hidden layers in the encoder for z.
                        (default: [512])
  --training-fraction TRAINING_FRACTION
                        Training detail: the fraction of the data used for
                        training. The rest is never seen by the inference
                        algorithm. Speeds up learning. (default: 0.9)
  --empty-drop-training-fraction FRACTION_EMPTIES
                        Training detail: the fraction of the training data
                        each epoch that is drawn (randomly sampled) from
                        surely empty droplets. (default: 0.2)
  --ignore-features BLACKLISTED_GENES [BLACKLISTED_GENES ...]
                        Integer indices of features to ignore entirely. In the
                        output count matrix, the counts for these features
                        will be unchanged. (default: [])
  --fpr FPR [FPR ...]   Target 'delta' false positive rate in [0, 1). Use 0
                        for a cohort of samples which will be jointly analyzed
                        for differential expression. A false positive is a
                        true signal count that is erroneously removed. More
                        background removal is accompanied by more signal
                        removal at high values of FPR. You can specify
                        multiple values, which will create multiple output
                        files. (default: [0.01])
  --exclude-feature-types EXCLUDE_FEATURES [EXCLUDE_FEATURES ...]
                        Feature types to ignore during the analysis. These
                        features will be left unchanged in the output file.
                        (default: [])
  --projected-ambient-count-threshold AMBIENT_COUNTS_IN_CELLS_LOW_LIMIT
                        Controls how many features are included in the
                        analysis, which can lead to a large speedup. If a
                        feature is expected to have less than
                        PROJECTED_AMBIENT_COUNT_THRESHOLD counts total in all
                        cells (summed), then that gene is excluded, and it
                        will be unchanged in the output count matrix. For
                        example, PROJECTED_AMBIENT_COUNT_THRESHOLD = 0 will
                        include all features which have even a single count in
                        any empty droplet. (default: 0.1)
  --learning-rate LEARNING_RATE
                        Training detail: lower learning rate for inference. A
                        OneCycle learning rate schedule is used, where the
                        upper learning rate is ten times this value. (For this
                        value, probably do not exceed 1e-3). (default: 0.0001)
  --checkpoint-mins CHECKPOINT_MIN
                        Checkpoint file will be saved periodically, with this
                        many minutes between each checkpoint. (default: 7.0)
  --final-elbo-fail-fraction FINAL_ELBO_FAIL_FRACTION
                        Training is considered to have failed if
                        (best_test_ELBO - final_test_ELBO)/(best_test_ELBO -
                        initial_test_ELBO) > FINAL_ELBO_FAIL_FRACTION.
                        Training will automatically re-run if --num-training-
                        tries > 1. By default, will not fail training based on
                        final_training_ELBO. (default: None)
  --epoch-elbo-fail-fraction EPOCH_ELBO_FAIL_FRACTION
                        Training is considered to have failed if
                        (previous_epoch_test_ELBO -
                        current_epoch_test_ELBO)/(previous_epoch_test_ELBO -
                        initial_train_ELBO) > EPOCH_ELBO_FAIL_FRACTION.
                        Training will automatically re-run if --num-training-
                        tries > 1. By default, will not fail training based on
                        epoch_training_ELBO. (default: None)
  --num-training-tries NUM_TRAINING_TRIES
                        Number of times to attempt to train the model. At each
                        subsequent attempt, the learning rate is multiplied by
                        LEARNING_RATE_RETRY_MULT. (default: 1)
  --learning-rate-retry-mult LEARNING_RATE_RETRY_MULT
                        Learning rate is multiplied by this amount each time a
                        new training attempt is made. (This parameter is only
                        used if training fails based on
                        EPOCH_ELBO_FAIL_FRACTION or FINAL_ELBO_FAIL_FRACTION
                        and NUM_TRAINING_TRIES is > 1.) (default: 0.2)
  --posterior-batch-size POSTERIOR_BATCH_SIZE
                        Training detail: size of batches when creating the
                        posterior. Reduce this to avoid running out of GPU
                        memory creating the posterior (will be slower).
                        (default: 128)
  --posterior-regularization {PRq,PRmu,PRmu_gene}
                        Posterior regularization method. (For experts: not
                        required for normal usage, see documentation). PRq is
                        approximate quantile-targeting. PRmu is approximate
                        mean-targeting aggregated over genes (behavior of
                        v0.2.0). PRmu_gene is approximate mean-targeting per
                        gene. (default: None)
  --alpha PRQ_ALPHA     Tunable parameter alpha for the PRq posterior
                        regularization method (not normally used: see
                        documentation). (default: None)
  --q CDF_THRESHOLD_Q   Tunable parameter q for the CDF threshold estimation
                        method (not normally used: see documentation).
                        (default: None)
  --estimator {map,mean,cdf,sample,mckp}
                        Output denoised count estimation method. (For experts:
                        not required for normal usage, see documentation).
                        (default: mckp)
  --estimator-multiple-cpu
                        Including the flag --estimator-multiple-cpu will use
                        more than one CPU to compute the MCKP output count
                        estimator in parallel (does nothing for other
                        estimators). (default: False)
  --constant-learning-rate
                        Including the flag --constant-learning-rate will use
                        the ClippedAdam optimizer instead of the OneCycleLR
                        learning rate schedule, which is the default. Learning
                        is faster with the OneCycleLR schedule. However,
                        training can easily be continued from a checkpoint for
                        more epochs than the initial command specified when
                        using ClippedAdam. On the other hand, if using the
                        OneCycleLR schedule with 150 epochs specified, it is
                        not possible to pick up from that final checkpoint and
                        continue training until 250 epochs. (default: False)
  --cpu-threads N_THREADS
                        Number of threads to use when pytorch is run on CPU.
                        Defaults to the number of logical cores. (default:
                        None)
  --debug               Including the flag --debug will log extra messages
                        useful for debugging. (default: False)
  --truth TRUTH_FILE    This is only used by developers for report generation.
                        Truth h5 file (for simulated data only). (default:
                        None)


```
