#　pySCENIC CLI

## STEP0 - パラメータSettings
```sh
 # 作業ディレクトリへ移動
  cd "$wd" || { echo "ERROR: Working directory '$wd' not found. Exiting." >&2; exit 1; }

  # Conda環境のアクティベート
  # この位置でアクティベートすることで、括弧内の全てのコマンドが指定された環境で実行されます。  
  # >>>>> ここに、あなたのConda初期化スクリプトのパスを追記 <<<<<
  # 例: source /home/tomohiro-nms/anaconda3/etc/profile.d/conda.sh
  source "$HOME/anaconda3/etc/profile.d/conda.sh" # <--- この行を追加・修正

  # Conda環境のアクティベート
  conda activate "$HOME/anaconda3/envs/pyscenic" # <--- この行はそのまま

  # --- 6. ジョブの開始メッセージと設定情報の出力 ---
  echo "--- pySCENIC 解析スクリプトを開始します ---"
  echo "開始日時: $(date)"
  echo "現在の作業ディレクトリ: $(pwd)"
  echo ""
  echo "--- ジョブ設定概要 ---"
  echo "入力発現行列 (Loom): ${INPUT_EXP_MTX}"
  echo "転写因子リスト (TF): ${TF_LIST}"
  echo "モチーフDBs: ${MOTIF_DB1}, ${MOTIF_DB2}" # ここで複数表示
  echo "アノテーションファイル: ${ANNOTATIONS_FNAME}"
  echo "出力ディレクトリ: ${OUTPUT_DIR}"
  echo "Daskワーカー数: ${NUM_WORKERS}"
  echo "Dask Memory Target: ${DASK_DISTRIBUTED__WORKER__MEMORY__TARGET:-N/A}"
  echo "Dask Memory Spill: ${DASK_DISTRIBUTED__WORKER__MEMORY__SPILL:-N/A}"
  echo "Dask Memory Pause: ${DASK_DISTRIBUTED__WORKER__MEMORY__PAUSE:-N/A}"
  echo "Dask Memory Terminate: ${DASK_DISTRIBUTED__WORKER__MEMORY__TERMINATE:-N/A}"
  echo "Dask Temporary Directory: ${DASK_TEMPORARY_DIRECTORY:-N/A}"
  echo "Python Warnings Filter: ${PYTHONWARNINGS:-N/A}"
  echo "-----------------------------------"
  echo ""

```


## STEP1 - GRN gene regulatory networkの構築
```sh
  # --- STEP 1: GRN構築 (Gene Regulatory Network Inference) ---
  echo -e "\n--- STEP 1/3: GRN構築を開始します (pyscenic grn) ---"

  if [ -s "${OUTPUT_ADJACENCIES_FILE}" ]; then
      echo "既存の ${OUTPUT_ADJACENCIES_FILE} が見つかりました。GRNステップをスキップします。"
  else
      echo "実行コマンド: pyscenic grn \"${INPUT_EXP_MTX}\" \"${TF_LIST}\" -o \"${OUTPUT_ADJACENCIES_FILE}\" --num_workers ${NUM_WORKERS} --method grnboost2 --transpose"
      pyscenic grn \
          "${INPUT_EXP_MTX}" \
          "${TF_LIST}" \
          -o "${OUTPUT_ADJACENCIES_FILE}" \
          --num_workers "${NUM_WORKERS}" \
          --method grnboost2

      # コマンドの終了ステータスを確認
      if [ ${PIPESTATUS[0]} -eq 0 ]; then
          echo ""
          echo "SUCCESS: pyscenic grn が正常に完了しました。"
      else
          echo "----------------------------------------------------"
          echo "ERROR: pyscenic grn の実行中にエラーが発生しました。" >&2
          echo "詳細については、上記の出力またはログファイルを確認してください。" >&2
          echo "終了ステータス: ${PIPESTATUS[0]}" >&2
          exit 1
      fi
  fi

```

## STEP2 - CTX
```sh
 # --- STEP 2: CTX構築 (Contextualization and motif enrichment) ---
  echo -e "\n--- STEP 2/3: CTX構築を開始します (pyscenic ctx) ---"

  if [ -s "${OUTPUT_REGULONS}" ]; then
      echo "既存の ${OUTPUT_REGULONS} が見つかりました。CTXステップをスキップします。"
  else
      echo "実行コマンド: pyscenic ctx \"${OUTPUT_ADJACENCIES_FILE}\" \"${MOTIF_DB1}\" \"${MOTIF_DB2}\" --annotations_fname \"${ANNOTATIONS_FNAME}\" --expression_mtx_fname \"${INPUT_EXP_MTX}\" --output \"${OUTPUT_REGULONS}\" --mask_dropouts --num_workers ${NUM_WORKERS} --mode \"dask_multiprocessing\""
      pyscenic ctx \
          "${OUTPUT_ADJACENCIES_FILE}" \
          "${MOTIF_DB1}" \
          "${MOTIF_DB2}" \
          --annotations_fname "${ANNOTATIONS_FNAME}" \
          --expression_mtx_fname "${INPUT_EXP_MTX}" \
          --output "${OUTPUT_REGULONS}" \
          --mask_dropouts \
          --num_workers "${NUM_WORKERS}" \
          --mode "dask_multiprocessing"

      if [ ${PIPESTATUS[0]} -eq 0 ]; then
          echo ""
          echo "SUCCESS: pyscenic ctx が正常に完了しました。"
      else
          echo "----------------------------------------------------"
          echo "ERROR: pyscenic ctx の実行中にエラーが発生しました。" >&2
          echo "詳細については、上記の出力またはログファイルを確認してください。" >&2
          echo "終了ステータス: ${PIPESTATUS[0]}" >&2
          exit 1
      fi
  fi
```

```sh
# --- STEP 3: AUCell (Activity calculation) ---
  echo -e "\n--- STEP 3/3: AUCellを開始します (pyscenic aucell) ---"

  if [ -s "${OUTPUT_AUC_LOOM}" ]; then
      echo "既存の ${OUTPUT_AUC_LOOM} が見つかりました。AUCellステップをスキップします。"
  else
      echo "実行コマンド: pyscenic aucell \"${INPUT_EXP_MTX}\" \"${OUTPUT_REGULONS}\" -o \"${OUTPUT_AUC_LOOM}\" --num_workers ${NUM_WORKERS}"
      pyscenic aucell \
          "${INPUT_EXP_MTX}" \
          "${OUTPUT_REGULONS}" \
          -o "${OUTPUT_AUC_LOOM}" \
          --num_workers "${NUM_WORKERS}"

      if [ ${PIPESTATUS[0]} -eq 0 ]; then
          echo ""
          echo "SUCCESS: pyscenic aucell が正常に完了しました。"
      else
          echo "----------------------------------------------------"
          echo "ERROR: pyscenic aucell の実行中にエラーが発生しました。" >&2
          echo "詳細については、上記の出力またはログファイルを確認してください。" >&2
          echo "終了ステータス: ${PIPESTATUS[0]}" >&2
          exit 1
      fi
  fi

  # --- ジョブの終了メッセージ ---
  echo ""
  echo "--- PySCENIC 解析スクリプトがすべて完了しました ---"
  echo "終了日時: $(date)"
  echo "-------------------------------------"


# 最終的なメッセージ（これはteeの範囲外なので、ターミナルに必ず表示されます）
echo "最終ログ: 全ての出力は '${LOG_FILE}' に保存されました。"

```



## STEP3 - AUCellの実行
