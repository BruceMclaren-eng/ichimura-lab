# リポジトリ整理案

## 基本方針

このリポジトリは、言語ごとではなく `プロジェクト単位` で整理する。

- 今の `practice/` は、1つの練習用プロジェクトとして独立させる
- 今後の研究は、別プロジェクトとして横に増やしていく
- GPU 性能解析の雛形や共通スクリプトは、各プロジェクトの外に共通資産として置く

## 整理後のトップレベル案

```text
ichimura-lab/
├─ README.md
├─ .gitignore
├─ plan.md
├─ projects/
│  ├─ practice-fem/
│  ├─ research-xxxx/
│  ├─ research-yyyy/
│  └─ shared/
├─ analysis/
│  └─ gpu/
│     ├─ formats/
│     ├─ scripts/
│     ├─ configs/
│     └─ templates/
├─ docs/
└─ scripts/
```

## 役割分担

### `projects/`

各研究・各練習を、独立した作業単位として置く場所。

ここには以下を入れる。

- そのプロジェクトのソースコード
- そのプロジェクト用のジョブスクリプト
- そのプロジェクトの実行結果
- そのプロジェクト固有のメモ
- そのプロジェクトの旧版や退避物

### `analysis/`

プロジェクト共通で使う解析資産を置く場所。

ここには以下を入れる。

- GPU 性能解析のフォーマット雛形
- roofline 用スクリプト
- profiler 用の共通設定
- CSV 整形や解析用の共通スクリプト
- 新しい研究でそのまま再利用したいテンプレート

つまり、

- プロジェクト固有の結果は `projects/<name>/`
- 使い回す解析資産は `analysis/`

で分ける。

## `practice/` の整理先

今の `practice/` は、まず `projects/practice-fem/` に移すのが自然。

想定構成は以下。

```text
projects/
└─ practice-fem/
   ├─ README.md
   ├─ src/
   │  ├─ fortran/
   │  ├─ cuda/
   │  └─ cpp/
   ├─ jobs/
   │  └─ fem_practice.sh
   ├─ results/
   │  ├─ profile/
   │  ├─ logs/
   │  └─ csv/
   ├─ archive/
   │  └─ dumped/
   └─ notes/
```

## `analysis/gpu/` の構成案

GPU 性能解析は今後も共通で使う前提なので、`projects/` と並列に置く。

```text
analysis/
└─ gpu/
   ├─ README.md
   ├─ formats/
   │  ├─ profiling_result_format.csv
   │  └─ gpu_efficiency_format.csv
   ├─ scripts/
   │  └─ roofline.py
   ├─ configs/
   │  ├─ a100_fp32.json
   │  └─ default_roofline_config.json
   └─ templates/
      └─ profiling-notes.md
```

### 各ディレクトリの意味

#### `formats/`

- 出力 CSV の雛形
- 列名の定義
- 共通の結果フォーマット

#### `scripts/`

- roofline 解析
- profiler 出力の整形
- 複数プロジェクトで再利用する解析コード

#### `configs/`

- GPU ごとの設定
- profiler 用の共通設定
- 実験条件テンプレート

#### `templates/`

- 新しい研究を始めるときの叩き台
- メモやレポートのテンプレート

## 今あるファイルの移し先

### `projects/practice-fem/src/`

- `practice/fem_solver/Fortran/*.f90`
- `practice/fem_solver/Fortran/matvec_cuda.cu`
- `practice/gputest/*.cu`
- `practice/gputest/*.f90`
- `practice/fem_solver/C++/practice.cpp`

### `projects/practice-fem/jobs/`

- `practice/job.sh`

### `projects/practice-fem/archive/`

- `practice/fem_solver/dumped/*`
- 比較用に残したい旧版ソース

### `projects/practice-fem/results/`

- 実行バイナリ
- `.exe`
- `.mod`
- `.o`
- `.ncu-rep`
- 実験で出たログや CSV

### `analysis/gpu/`

以下は、共通化できるものを `analysis/gpu/` に寄せる。

- `practice/fem_solver/Gpu_analysis/analysis_format/roofline.py`
- `practice/fem_solver/Gpu_analysis/analysis_format/config.json`

以下は中身を見て判断する。

- `practice/fem_solver/Gpu_analysis/Projects/2dlaplace_fem/roofline.py`
  - 共通化できるなら `analysis/gpu/scripts/`
  - その実験専用なら `projects/practice-fem/notes/` か `results/`
- `practice/fem_solver/Gpu_analysis/Projects/2dlaplace_fem/config.json`
  - 汎用設定なら `analysis/gpu/configs/`
  - 実験専用ならプロジェクト側に残す
- `practice/fem_solver/Gpu_analysis/output/*.csv`
  - 雛形なら `analysis/gpu/formats/`
  - 実結果なら `projects/practice-fem/results/csv/`

## `.gitignore` の方針

ルートに `.gitignore` を置いて、少なくとも以下を除外する。

```gitignore
# local settings
.claude/

# per-project outputs
projects/*/results/

# compiled artifacts
*.exe
*.out
*.mod
*.o
*.obj
*.a
*.so
*.dll
*.exp
*.lib

# profiler outputs
*.ncu-rep

# logs
*.log
```

必要なら、手で整えた軽量な結果だけ個別に追跡する。

## 移行の進め方

### Phase 1

1. ルートに `README.md` を追加
2. `projects/` を作成
3. `analysis/gpu/` を作成
4. `.gitignore` を追加

### Phase 2

1. `projects/practice-fem/` を作成
2. `src/`, `jobs/`, `results/`, `archive/`, `notes/` を作成
3. `practice/` の中身を役割ごとに移動

### Phase 3

1. `Gpu_analysis/` から共通資産を `analysis/gpu/` に切り出す
2. 実験結果だけを `projects/practice-fem/results/` に残す
3. `analysis/gpu/README.md` を追加して、再利用方法を書く

### Phase 4

1. 今後の研究を `projects/research-.../` として追加
2. GPU 解析は `analysis/gpu/` を共通利用
3. 本当に共通化が必要なものだけ `projects/shared/` に出す

## この案の利点

- `practice` を練習用プロジェクトとしてきれいに隔離できる
- 今後の研究を追加しやすい
- GPU 解析の雛形をすぐ再利用できる
- プロジェクト固有の成果物と共通資産が混ざらない

## 次にやること

この方針でよければ、次は以下を実際に進める。

1. `README.md` の追加
2. `.gitignore` の追加
3. `projects/practice-fem/` の雛形作成
4. `analysis/gpu/` の雛形作成
5. `practice/` と `Gpu_analysis/` の移設
