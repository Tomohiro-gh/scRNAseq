#　pySCENIC CLI







## STEP1 - GRN gene regulatory networkの構築


## STEP2 - CTX
```sh
#!/bin/bash

# --- 設定変数 ---
# ご自身の環境に合わせてパスを修正してください
COEXP_MATRIX="/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp221_5tissues_DEG_Capillaries/X_youngCaps_pySCENIC/results_pyscenic_small/adjacencies.tsv"
ANNOTATIONS_FNAME="//Users/tomohiro/Library/CloudStorage/Dropbox/Database/SCENIC/Auxiliary_Dataset_mm10_mc_v10_clust/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
FEATHER_FNAME="/Users/tomohiro/Library/CloudStorage/Dropbox/Database/SCENIC/Auxiliary_Dataset_mm10_mc_v10_clust/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
EXPRESSION_MTX_FNAME="/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp221_5tissues_DEG_Capillaries/X_youngCaps_pySCENIC/young_cap_small.loom"
OUTPUT_REGULONS="/Users/tomohiro/Library/CloudStorage/Dropbox/FukuharaLab_Res/Experiment/Exp221_5tissues_DEG_Capillaries/X_youngCaps_pySCENIC/results_pyscenic_small/regulons.csv"
NUM_WORKERS=12


# --- pyscenic ctx コマンドの実行 ---
echo ""
echo "pyscenic ctx の実行を開始します..."
echo "コマンド: pyscenic ctx ${COEXP_MATRIX} --annotations_fname ${ANNOTATIONS_FNAME} --expression_mtx_fname ${EXPRESSION_MTX_FNAME} --output ${OUTPUT_REGULONS} --num_workers ${NUM_WORKERS} --mode \"dask_multiprocessing\""
echo "----------------------------------------------------"

# コマンドを実行し、標準エラー出力 (2) を標準出力 (1) にリダイレクト
# tee を使って、画面にも出力しつつ、ログファイルにも保存する
pyscenic ctx \
    "${COEXP_MATRIX}" \
    "${FEATHER_FNAME}" \
    --annotations_fname "${ANNOTATIONS_FNAME}" \
    --expression_mtx_fname "${EXPRESSION_MTX_FNAME}" \
    --output "${OUTPUT_REGULONS}" \
    --num_workers "${NUM_WORKERS}" \
    --mode "dask_multiprocessing" 2>&1 | tee pyscenic_ctx_log.txt

# コマンドの終了ステータスを確認
if [ ${PIPESTATUS[0]} -eq 0 ]; then
    echo "----------------------------------------------------"
    echo "SUCCESS: pyscenic ctx が正常に完了しました。"
else
    echo "----------------------------------------------------"
    echo "ERROR: pyscenic ctx の実行中にエラーが発生しました。" >&2
    echo "詳細については、上記の出力または 'pyscenic_ctx_log.txt' を確認してください。" >&2
    exit 1
fi
```





## STEP3 - AUCellの実行
