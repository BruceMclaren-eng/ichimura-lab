# practice-fem

FEM Solver と GPU/OpenACC 関連の練習用プロジェクト。

## 構成

- `src/`: ソースコード
- `jobs/`: 実行ジョブ
- `results/`: 実行結果、バイナリ、プロファイル出力
- `archive/`: 退避した旧版や比較用ファイル
- `notes/`: 実験メモ

## 注意

- `results/` 配下は原則 Git で追跡しない
- 共通の GPU 解析スクリプトは `../../analysis/gpu/` を使う
