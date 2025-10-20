#!/usr/bin/env python3
"""
merge_docking_results.py

功能说明：
    本脚本用于将 MEGADOCK、HDOCK、Alphafold 的对接结果合并为一个统一的表格。
    输出文件包含五列：
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
"""

import pandas as pd
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Merge MEGADOCK, HDOCK, and AlphaFold results into one table."
    )
    parser.add_argument("--megadock", required=True, help="Path to megadock.tsv file")
    parser.add_argument("--hdock", required=True, help="Path to hdock.tsv file")
    parser.add_argument("--af", required=True, help="Path to AlphaFold af_c.tsv file")
    parser.add_argument("--output", required=True, help="Output TSV file name")

    args = parser.parse_args()

    try:
        # 读取 MEGADOCK 文件
        megadock = pd.read_csv(args.megadock, sep="\t", header=None, names=["ID1", "ID2", "MEGADOCK_Score"])

        # 读取 HDOCK 文件
        hdock = pd.read_csv(args.hdock, sep="\t", header=None, names=["ID1", "ID2", "HDOCK_Score"])

        # 读取 AlphaFold 文件
        af = pd.read_csv(args.af, sep="\t")
        af[['ID1', 'ID2']] = af['Pair'].str.split('-', expand=True)
        af['Alphafold_pTM+ipTM'] = af['PTM'] + af['IPTM']
        af = af[['ID1', 'ID2', 'Alphafold_pTM+ipTM']]

        # 合并三个文件
        df = pd.merge(megadock, hdock, on=['ID1', 'ID2'], how='outer')
        df = pd.merge(df, af, on=['ID1', 'ID2'], how='outer')

        # 输出结果
        df.to_csv(args.output, sep="\t", index=False)
        print(f"✅ 合并完成，输出文件：{args.output}")

    except Exception as e:
        print(f"❌ 出错：{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
