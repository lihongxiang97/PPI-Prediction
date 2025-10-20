#!/usr/bin/env python3
"""
merge_docking_results.py

功能说明：
    本脚本用于将 MEGADOCK、HDOCK、Alphafold 的对接结果合并为一个统一的表格。
    输出文件包含以下列（视输入文件而定）：
        ID1    ID2    MEGADOCK_Score    HDOCK_Score    Alphafold_pTM+ipTM

输入文件格式：
    1. megadock.tsv:  (无表头)
        ID1    ID2    MEGADOCK_Score
    2. hdock.tsv:     (无表头)
        ID1    ID2    HDOCK_Score
    3. af_c.tsv:      (有表头)
        Pair    PTM    IPTM
        （Pair形如 AT1G01010-AT1G01110）

使用方法：
    python merge_docking_results.py \
        --megadock megadock.tsv \
        --hdock hdock.tsv \
        --af af_c.tsv \
        --output merged_results.tsv

说明：
    若只提供一个或两个输入文件，脚本会自动合并可用的部分并输出。
"""

import pandas as pd
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(
        description="Merge MEGADOCK, HDOCK, and AlphaFold results into one table."
    )
    parser.add_argument("--megadock", help="Path to megadock.tsv file (optional)")
    parser.add_argument("--hdock", help="Path to hdock.tsv file (optional)")
    parser.add_argument("--af", help="Path to AlphaFold af_c.tsv file (optional)")
    parser.add_argument("--output", required=True, help="Output TSV file name")

    args = parser.parse_args()

    # 检查至少有一个输入文件存在
    if not any([args.megadock, args.hdock, args.af]):
        print("❌ 错误：请至少提供一个输入文件 (--megadock / --hdock / --af)")
        sys.exit(1)

    dfs = []  # 存放可用的 DataFrame
    merge_cols = ['ID1', 'ID2']

    try:
        if args.megadock and os.path.exists(args.megadock):
            print(f"📄 读取 MEGADOCK 文件: {args.megadock}")
            megadock = pd.read_csv(args.megadock, sep="\t", header=None, names=["ID1", "ID2", "MEGADOCK_Score"])
            dfs.append(megadock)
        if args.hdock and os.path.exists(args.hdock):
            print(f"📄 读取 HDOCK 文件: {args.hdock}")
            hdock = pd.read_csv(args.hdock, sep="\t", header=None, names=["ID1", "ID2", "HDOCK_Score"])
            dfs.append(hdock)
        if args.af and os.path.exists(args.af):
            print(f"📄 读取 AlphaFold 文件: {args.af}")
            af = pd.read_csv(args.af, sep="\t")
            af[['ID1', 'ID2']] = af['Pair'].str.split('-', expand=True)
            af['Alphafold_pTM+ipTM'] = (af['PTM'] + af['IPTM']).round(2)
            af = af[['ID1', 'ID2', 'Alphafold_pTM+ipTM']]
            dfs.append(af)

        # 合并所有可用的文件
        if len(dfs) == 1:
            merged = dfs[0]
        else:
            merged = dfs[0]
            for df in dfs[1:]:
                merged = pd.merge(merged, df, on=merge_cols, how='outer')

        # 按 ID 排序（可选）
        merged = merged.sort_values(by=["ID1", "ID2"], ignore_index=True)

        # 保存结果
        merged.to_csv(args.output, sep="\t", index=False)
        print(f"✅ 合并完成，输出文件：{args.output}")

    except Exception as e:
        print(f"❌ 出错：{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
