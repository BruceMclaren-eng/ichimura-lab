# Roofline Analysis Tool
# ncu_report ライブラリを使って .ncu-rep からカーネルごとのRoofline指標をCSV出力する
#
# 使い方:
#   set PYTHONPATH=C:\Program Files\NVIDIA Corporation\Nsight Compute 2026.1.1\extras\python
#   python roofline.py --rep roop1_profile..ncu-rep --config config.json --output output/roofline.csv

import argparse
import json
import os
import warnings
import pandas as pd
import ncu_report


def get_metric(kernel, metric_name):
    """メトリクス値を取得。存在しない場合は 0.0 を返す。"""
    try:
        return kernel[metric_name].value()
    except Exception:
        return 0.0


def find_config_note(kernel_name, kernel_cfg):
    """
    config の kernels キーと ncu-rep の Mangled Name を部分一致で照合し、
    一致したエントリの note を返す。一致しない場合は空文字。
    """
    for cfg_key, cfg_val in kernel_cfg.items():
        if cfg_key in kernel_name or kernel_name in cfg_key:
            return cfg_val.get("note", "")
    return ""


def analyze(rep_path, config_path, output_path):

    # --- config 読み込み ---
    with open(config_path, "r", encoding="utf-8") as f:
        config = json.load(f)

    gpu_cfg    = config.get("gpu", {})
    kernel_cfg = config.get("kernels", {})

    precision        = gpu_cfg.get("precision", "fp64")
    peak_fp32_gflops = gpu_cfg.get("peak_fp32_gflops", 19500)   # GFLOPS
    peak_fp64_gflops = gpu_cfg.get("peak_fp64_gflops", 9700)    # GFLOPS
    peak_dram_bw     = gpu_cfg.get("peak_dram_bandwidth_gbs", 1555)  # GB/s

    peak_flops  = peak_fp64_gflops if precision == "fp64" else peak_fp32_gflops
    ridge_point = peak_flops / peak_dram_bw  # FLOP/Byte

    # --- ncu-rep 読み込み ---
    print(f"[INFO] レポートを読み込み中: {rep_path}")
    report = ncu_report.load_report(rep_path)

    results = []

    for range_idx in range(report.num_ranges()):
        r = report[range_idx]
        for kernel_idx in range(len(r)):
            kernel = r[kernel_idx]

            # Function Name: ncu-rep のカーネル名（Mangled Name）
            function_name = kernel.name()

            # Kernel Name: config の note（Mangled Name でキー照合）
            kernel_note = find_config_note(function_name, kernel_cfg)

            # --- precisionに応じたメトリクス名を選択 ---
            cycles_per_sec = get_metric(kernel, "smsp__cycles_elapsed.avg.per_second")

            if precision == "fp32":
                fadd_per_cycle = get_metric(
                    kernel, "smsp__sass_thread_inst_executed_op_fadd_pred_on.sum.per_cycle_elapsed"
                )
                fmul_per_cycle = get_metric(
                    kernel, "smsp__sass_thread_inst_executed_op_fmul_pred_on.sum.per_cycle_elapsed"
                )
                ffma_x2 = get_metric(
                    kernel, "derived__smsp__sass_thread_inst_executed_op_ffma_pred_on_x2"
                )
            else:  # fp64
                fadd_per_cycle = get_metric(
                    kernel, "smsp__sass_thread_inst_executed_op_dadd_pred_on.sum.per_cycle_elapsed"
                )
                fmul_per_cycle = get_metric(
                    kernel, "smsp__sass_thread_inst_executed_op_dmul_pred_on.sum.per_cycle_elapsed"
                )
                ffma_x2 = get_metric(
                    kernel, "derived__smsp__sass_thread_inst_executed_op_dfma_pred_on_x2"
                )

            # --- 各 Work [GFLOPS] ---
            fadd_gflops  = fadd_per_cycle * cycles_per_sec * 1e-9
            fmul_gflops  = fmul_per_cycle * cycles_per_sec * 1e-9
            ffma_gflops  = ffma_x2        * cycles_per_sec * 1e-9
            total_gflops = fadd_gflops + fmul_gflops + ffma_gflops

            # --- Achieved Traffic [GBytes/sec] ---
            achieved_traffic = get_metric(kernel, "dram__bytes.sum.per_second") * 1e-9

            # --- Arithmetic Intensity [FLOP/Byte] ---
            ai = total_gflops / achieved_traffic if achieved_traffic > 0 else 0.0

            # --- Bottleneck Type ---
            bottleneck = "Memory Bound" if ai < ridge_point else "Compute Bound"

            # --- Theoretical Peak Performance [GFLOPS] ---
            if bottleneck == "Memory Bound":
                theoretical_peak = peak_dram_bw * ai
            else:
                theoretical_peak = peak_flops

            # --- Achieved Efficiency [%] ---
            achieved_efficiency = (
                (total_gflops / theoretical_peak * 100) if theoretical_peak > 0 else 0.0
            )

            results.append({
                "Function Name":                      function_name,
                "Kernel Name":                        kernel_note,
                "Achieved Fadd Work [GFLOPS]":        fadd_gflops,
                "Achieved Fmul Work [GFLOPS]":        fmul_gflops,
                "Achieved Ffma Work [GFLOPS]":        ffma_gflops,
                "Total Achieved Work [GFLOPS]":       total_gflops,
                "Achieved Traffic [GBytes/sec]":      achieved_traffic,
                "Arithmetic Intensity [FLOP/Bytes]":  ai,
                "Ridge Point [FLOP/Bytes]":           ridge_point,
                "Bottleneck Type":                    bottleneck,
                "Theoretical Peak Performance [GFLOPS]": theoretical_peak,
                "Achieved Efficiency [%]":            achieved_efficiency,
            })

    # --- DataFrame 化・CSV 出力 ---
    df = pd.DataFrame(results)

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    df.to_csv(output_path, index=False, float_format="%.6f")
    print(f"[INFO] 結果を保存しました: {output_path}")

    # --- コンソール表示 ---
    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 300)
    pd.set_option("display.float_format", "{:.4f}".format)
    print(df.to_string(index=False))

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Roofline analysis from ncu-rep")
    parser.add_argument("--rep",    required=True,                   help=".ncu-rep ファイルのパス")
    parser.add_argument("--config", required=True,                   help="config.json のパス")
    parser.add_argument("--output", default="output/roofline.csv",   help="出力CSVのパス")
    args = parser.parse_args()

    analyze(args.rep, args.config, args.output)
