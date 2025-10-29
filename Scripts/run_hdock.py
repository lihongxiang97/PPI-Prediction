#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Hongxiang Li, 2025/07/09 (Modified: add pdbs symlink + parallel support + lowercase pair list)

import os
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed


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

    # 判断长度，决定谁是 R/L
    len1 = get_pdb_length(pdb1)
    len2 = get_pdb_length(pdb2)
    R, L = (id1, id2) if len1 >= len2 else (id2, id1)

    # 输出文件名
    out_name = os.path.join(output_dir, f"{R}-{L}.out")
    out_pdb = os.path.join(output_dir, f"{R}-{L}.out.pdb")

    # 准备命令路径
    hdock_cmd = f"{hdock_path}/hdock" if hdock_path else "hdock"
    createpl_cmd = f"{hdock_path}/createpl" if hdock_path else "createpl"

    cwd = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)

    # 创建软链接到 pdb_dir
    symlink_path = os.path.join(output_dir, "pdbs")
    if not os.path.exists(symlink_path):
        try:
            os.symlink(pdb_dir, symlink_path)
            print(f"[INFO] Linked pdbs -> {pdb_dir}")
        except FileExistsError:
            pass

    os.chdir(output_dir)

    try:
        # === 跳过逻辑 ===
        if os.path.exists(out_pdb):
            print(f"[SKIP] {R}-{L}: Final PDB exists, skip all.")
            # 直接提取得分
            with open(out_pdb) as f:
                for line in f:
                    if line.startswith("REMARK Score:"):
                        score = line.strip().replace("REMARK Score: ", "")
                        return R, L, float(score)
            return None

        elif os.path.exists(out_name):
            print(f"[RESUME] {R}-{L}: Out file exists, run createpl only.")
            subprocess.run([
                createpl_cmd,
                f"{R}-{L}.out",
                f"{R}-{L}.out.pdb",
                "-nmax", "1",
                "-complex"
            ], check=True)

        else:
            print(f"[RUN] {R}-{L}: Running hdock + createpl.")
            subprocess.run([
                hdock_cmd,
                f"pdbs/{R}.pdb",
                f"pdbs/{L}.pdb",
                "-spacing", "1.2",
                "-angle", "15",
                "-out", f"{R}-{L}.out"
            ], check=True, text=True)

            subprocess.run([
                createpl_cmd,
                f"{R}-{L}.out",
                f"{R}-{L}.out.pdb",
                "-nmax", "1",
                "-complex"
            ], check=True)

        # 提取得分
        if os.path.exists(f"{R}-{L}.out.pdb"):
            with open(f"{R}-{L}.out.pdb") as f:
                for line in f:
                    if line.startswith("REMARK Score:"):
                        score = line.strip().replace("REMARK Score: ", "")
                        return R, L, float(score)

    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed for {R}-{L}: {e}")
    except Exception as e:
        print(f"[ERROR] Unexpected error for {R}-{L}: {e}")
    finally:
        os.chdir(cwd)

    return None


def main():
    parser = argparse.ArgumentParser(description="Run HDOCK for protein pairs in parallel")
    parser.add_argument("-l","--pair_list", required=True, help="Protein pair list file")
    parser.add_argument("-d","--pdb_dir", required=True, help="Directory with PDB files")
    parser.add_argument("-od","--output_dir", required=True, help="Directory to store HDOCK outputs")
    parser.add_argument("-p","--hdock_path", default=None, help="Optional path to HDOCK executables")
    parser.add_argument("-r","--result_file", default="hdock_result.txt", help="Output result file")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of parallel HDOCK tasks (default: 8)")
    args = parser.parse_args()

    args.output_dir = os.path.abspath(args.output_dir)
    args.pdb_dir = os.path.abspath(args.pdb_dir)
    args.result_file = os.path.abspath(args.result_file)
    args.hdock_path = os.path.abspath(args.hdock_path) if args.hdock_path else None
    
    # ✅ 读取所有配对，并将字母转为小写
    pairs = []
    with open(args.pair_list, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("ID"):
                id1, id2 = line.strip().split()
                id1, id2 = id1.lower(), id2.lower()   # ✅ 转为小写
                pairs.append((id1, id2))

    print(f"[INFO] Loaded {len(pairs)} pairs. Running with {args.threads} threads...")

    results = []

    # 并行运行 HDOCK
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        future_to_pair = {
            executor.submit(run_hdock_on_pair, id1, id2, args.pdb_dir, args.output_dir, args.hdock_path): (id1, id2)
            for id1, id2 in pairs
        }

        for future in as_completed(future_to_pair):
            pair = future_to_pair[future]
            try:
                result = future.result()
                if result:
                    results.append(result)
                    R, L, score = result
                    print(f"[OK] {R}-{L}: Score={score:.3f}")
            except Exception as e:
                print(f"[ERROR] Exception in {pair}: {e}")

    # 排序并写出结果
    results.sort(key=lambda x: x[2], reverse=True)
    with open(args.result_file, "w") as out:
        for R, L, score in results:
            out.write(f"{R}\t{L}\t{score:.4f}\n")

    print(f"[DONE] All {len(pairs)} pairs processed. Results saved to {args.result_file}")


if __name__ == "__main__":
    main()
