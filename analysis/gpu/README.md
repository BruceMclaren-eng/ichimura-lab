# GPU Analysis

複数プロジェクトで使い回す GPU 性能解析資産を置く場所。

## 構成

- `formats/`: CSV 形式や出力列の雛形
- `scripts/`: roofline 解析などの共通スクリプト
- `configs/`: GPU や profiler の共通設定
- `templates/`: 新規プロジェクト向けテンプレート

## 使い方

- 各プロジェクトの profiler 出力は、そのプロジェクトの `results/` に置く
- 解析スクリプトや共通 config はここから参照する
- プロジェクト固有の解析結果までここに置かない
