#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Hongxiang Li, 2025/07/09 (modified: added .out existence check)

import os
import argparse
import subprocess
import sys

def get_pdb_length(pdb_file):
    """计算PDB中氨基酸的数量（CA原子行数）"""
    count = 0
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                count += 1
    return count

def run_megadock(pdb1, pdb2, output_dir, pdb_dir, docker_image, n_decoys, t, cpu_cores):
    """执行 MEGADOCK 对接和得分计算"""
    id1 = os.path.basename(pdb1).replace(".pdb", "")
    id2 = os.path.basename(pdb2).replace(".pdb", "")

    # 判断长短，决定 -R -L
    len1 = get_pdb_length(pdb1)
    len2 = get_pdb_length(pdb2)

    if len1 >= len2:
        R, L = id1, id2
    else:
        R, L = id2, id1

    # 准备文件路径
    os.makedirs(output_dir, exist_ok=True)
    out_basename = f"{R}-{L}.out"
    out_path_host = os.path.abspath(os.path.join(output_dir, out_basename))

    # 如果 .out 已存在，跳过对接
    if os.path.exists(out_path_host):
        print(f"[INFO] Output {out_basename} already exists, skipping docking step.")
    else:
        # MEGADOCK 对接命令
        dock_cmd = (
            "docker run --rm --gpus all "
            f"-e OMP_NUM_THREADS={cpu_cores} "
            f"-v {os.path.abspath(pdb_dir)}:/opt/MEGADOCK/data "
            f"-v {os.path.abspath(output_dir)}:/opt/MEGADOCK/out "
            f"{docker_image} "
            "megadock-gpu "
            f"-R /opt/MEGADOCK/data/{R}.pdb "
            f"-L /opt/MEGADOCK/data/{L}.pdb "
            f"-o /opt/MEGADOCK/out/{out_basename} "
            f"-N {n_decoys} -t {t}"
        )
        print(f"[INFO] Running MEGADOCK for {R} vs {L}") 
        ret = subprocess.run(dock_cmd, shell=True)
        if ret.returncode != 0:
            print(f"[ERROR] MEGADOCK docking failed for {R} vs {L} (return code {ret.returncode})")
            return R, L, None

    # 获取得分
    ppiscore_cmd = (
        "docker run --rm "
        f"-v {os.path.abspath(output_dir)}:/opt/MEGADOCK/out "
        f"{docker_image} "
        f"ppiscore out/{out_basename} {n_decoys}"
    )
    try:
        result = subprocess.check_output(ppiscore_cmd, shell=True, stderr=subprocess.STDOUT).decode()
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] ppiscore failed for {R} vs {L}: {e.output.decode().strip()}")
        return R, L, None

    # 提取得分
    for line in result.splitlines():
        if ".out" in line and "E =" in line:
            # 例: something, E = -123.45, ...
            parts = line.split(",")
            for part in parts:
                if "E =" in part:
                    score_str = part.strip().replace("E =", "").strip()
                    try:
                        return R, L, float(score_str)
                    except ValueError:
                        print(f"[WARNING] Unable to parse score '{score_str}' for {R} vs {L}")
                        return R, L, None

    print(f"[WARNING] No score parsed for {R} vs {L} (ppiscore output did not contain expected pattern)")
    return R, L, None

def main():
    parser = argparse.ArgumentParser(description="Run MEGADOCK for PPI prediction")
    parser.add_argument("--pair_list", required=True, help="Protein pair list file (ID1 ID2)")
    parser.add_argument("--pdb_dir", required=True, help="Directory containing PDB files")
    parser.add_argument("--output_dir", required=True, help="MEGADOCK output directory")
    parser.add_argument("--result_file", default="megadock_result.txt", help="File to save MEGADOCK scores")
    parser.add_argument("--docker_image", default="hub.rat.dev/akiyamalab/megadock:gpu", help="Docker image for MEGADOCK")
    parser.add_argument("-N", type=int, default=10800, help="Number of decoys, default 10800")
    parser.add_argument("-t", type=int, default=3, help="Thread number for FFT, default 3")
    parser.add_argument("-e", type=int, default=32, help="Number of CPU cores (OMP_NUM_THREADS), default 32")

    args = parser.parse_args()

    results = []

    if not os.path.isfile(args.pair_list):
        print(f"[ERROR] Pair list file not found: {args.pair_list}")
        sys.exit(1)

    with open(args.pair_list, 'r') as f:
        for line in f:
            if line.strip() and not line.strip().startswith("ID"):
                parts = line.strip().split()
                if len(parts) < 2:
                    print(f"[WARNING] Skipping malformed line: {line.strip()}")
                    continue
                id1, id2 = parts[0], parts[1]
                pdb1 = os.path.join(args.pdb_dir, f"{id1}.pdb")
                pdb2 = os.path.join(args.pdb_dir, f"{id2}.pdb")

                if not os.path.exists(pdb1) or not os.path.exists(pdb2):
                    print(f"[WARNING] Missing PDB file for pair: {id1}, {id2}")
                    continue

                R, L, score = run_megadock(
                    pdb1, pdb2,
                    output_dir=args.output_dir,
                    pdb_dir=args.pdb_dir,
                    docker_image=args.docker_image,
                    n_decoys=args.N,
                    t=args.t,
                    cpu_cores=args.e
                )

                if score is not None:
                    results.append(f"{R}\t{L}\t{score:.4f}")

    # 按得分从大到小排序
    results_sorted = sorted(results, key=lambda x: float(x.strip().split("\t")[2]), reverse=True)

    # 写入输出文件
    try:
        with open(args.result_file, 'w') as out:
            out.write("\n".join(results_sorted))
        print(f"[DONE] MEGADOCK completed. Results saved to {args.result_file}")
    except Exception as e:
        print(f"[ERROR] Failed to write result file {args.result_file}: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
