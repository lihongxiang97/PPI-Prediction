#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Hongxiang Li, 2025/07/09

import os
import argparse
import subprocess

def get_pdb_length(pdb_file):
    """统计PDB中CA原子的数量"""
    count = 0
    with open(pdb_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") and " CA " in line:
                count += 1
    return count

def run_hdock_on_pair(id1, id2, pdb_dir, output_dir, hdock_path=None):
    pdb1 = os.path.join(pdb_dir, f"{id1}.pdb")
    pdb2 = os.path.join(pdb_dir, f"{id2}.pdb")

    if not os.path.exists(pdb1) or not os.path.exists(pdb2):
        print(f"[WARNING] Missing PDB file for {id1} or {id2}")
        return None

    # 判断长度
    len1 = get_pdb_length(pdb1)
    len2 = get_pdb_length(pdb2)

    if len1 >= len2:
        R, L = id1, id2
    else:
        R, L = id2, id1

    # 输出文件名（保持ID顺序，用户可识别）
    out_name = f"{R}-{L}.out"
    out_pdb = f"{out_name}.pdb"

    # HDOCK 命令
    hdock_cmd = f"{hdock_path}/hdock" if hdock_path else "hdock"
    createpl_cmd = f"{hdock_path}/createpl" if hdock_path else "createpl"

    # 切换目录执行 HDOCK
    cwd = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)

    try:
        subprocess.run([
            hdock_cmd,
            f"{pdb_dir}/{R}.pdb",
            f"{pdb_dir}/{L}.pdb",
            "-spacing", "1.2",
            "-angle", "15",
            "-out", out_name
        ], check=True)

        subprocess.run([
            createpl_cmd,
            out_name,
            out_pdb,
            "-nmax", "1",
            "-complex"
        ], check=True)

        # 提取得分
        if os.path.exists(out_pdb):
            with open(out_pdb) as f:
                for line in f:
                    if line.startswith("REMARK Score:"):
                        score = line.strip().replace("REMARK Score: ", "")
                        return R, L, float(score)
    except Exception as e:
        print(f"[ERROR] Failed for {id1}-{id2}: {e}")
    finally:
        os.chdir(cwd)

    return None

def main():
    parser = argparse.ArgumentParser(description="Run HDOCK for protein pairs")
    parser.add_argument("--pair_list", required=True, help="Protein pair list file")
    parser.add_argument("--pdb_dir", required=True, help="Directory with PDB files")
    parser.add_argument("--output_dir", required=True, help="Directory to store HDOCK outputs")
    parser.add_argument("--hdock_path", default=None, help="Optional path to HDOCK executables")
    parser.add_argument("--result_file", default="hdock_result.txt", help="Output result file")

    args = parser.parse_args()

    results = []

    with open(args.pair_list, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("ID"):
                id1, id2 = line.strip().split()
                result = run_hdock_on_pair(id1, id2, args.pdb_dir, args.output_dir, args.hdock_path)
                if result:
                    R, L, score = result
                    results.append((R, L, score))

    # 按得分从高到低排序
    results.sort(key=lambda x: x[2], reverse=True)

    # 写入结果
    with open(args.result_file, "w") as out:
        for R, L, score in results:
            out.write(f"{R}\t{L}\t{score:.4f}\n")

    print(f"[DONE] HDOCK completed. Results saved to {args.result_file}")

if __name__ == "__main__":
    main()
