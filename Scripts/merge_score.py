#!/usr/bin/env python3
"""
merge_docking_results.py

功能说明：
    将 MEGADOCK、HDOCK、AlphaFold 的结果文件合并为统一表格。
    输出列（根据输入文件而定）：
        ID1    ID2    MEGADOCK_Score    HDOCK_Score    Alphafold_pTM+ipTM

特性：
    ✅ 自动将 ID 转为大写
    ✅ 自动忽略 ID 顺序（AT1G01010-AT1G01110 == AT1G01110-AT1G01010）
    ✅ 可仅输入 1、2 或 3 个文件
    ✅ 结果自动按 ID 排序并输出 TSV

输入文件格式：
    1. megadock.tsv: 无表头，三列：ID1 ID2 MEGADOCK_Score
    2. hdock.tsv:    无表头，三列：ID1 ID2 HDOCK_Score
    3. af_c.tsv:     有表头：Pair PTM IPTM （Pair 形如 AT1G01010-AT1G01110）

使用示例：
    python merge_docking_results.py \
        --megadock megadock.tsv \
        --hdock hdock.tsv \
        --af af_c.tsv \
        --output merged_results.tsv
"""

import pandas as pd
import argparse
import sys
import os


def normalize_ids(df, id1_col="ID1", id2_col="ID2"):
    """将ID转为大写并按字母顺序重新排列，保证ID1 < ID2"""
    df[id1_col] = df[id1_col].str.upper()
    df[id2_col] = df[id2_col].str.upper()
    df[['ID1', 'ID2']] = pd.DataFrame(
        df[[id1_col, id2_col]].apply(lambda x: sorted([x[id1_col], x[id2_col]]), axis=1).tolist(),
        index=df.index
    )
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Merge MEGADOCK, HDOCK, and AlphaFold results into one table."
    )
    parser.add_argument("--megadock", help="Path to megadock.tsv file (optional)")
    parser.add_argument("--hdock", help="Path to hdock.tsv file (optional)")
    parser.add_argument("--af", help="Path to AlphaFold af_c.tsv file (optional)")
    parser.add_argument("--output", required=True, help="Output TSV file name")

    args = parser.parse_args()

    if not any([args.megadock, args.hdock, args.af]):
        print("❌ 错误：请至少提供一个输入文件 (--megadock / --hdock / --af)")
        sys.exit(1)

    dfs = []
    merge_cols = ['ID1', 'ID2']

    try:
        if args.megadock and os.path.exists(args.megadock):
            print(f"📄 读取 MEGADOCK 文件: {args.megadock}")
            megadock = pd.read_csv(args.megadock, sep="\t", header=None, names=["ID1", "ID2", "MEGADOCK_Score"])
            megadock = normalize_ids(megadock)
            dfs.append(megadock)

        if args.hdock and os.path.exists(args.hdock):
            print(f"📄 读取 HDOCK 文件: {args.hdock}")
            hdock = pd.read_csv(args.hdock, sep="\t", header=None, names=["ID1", "ID2", "HDOCK_Score"])
            hdock = normalize_ids(hdock)
            dfs.append(hdock)

        if args.af and os.path.exists(args.af):
            print(f"📄 读取 AlphaFold 文件: {args.af}")
            af = pd.read_csv(args.af, sep="\t")
            af[['ID1', 'ID2']] = af['Pair'].str.split('-', expand=True)
            af = normalize_ids(af)
            af['Alphafold_pTM+ipTM'] = (af['PTM'] + af['IPTM']).round(2)
            af = af[['ID1', 'ID2', 'Alphafold_pTM+ipTM']]
            dfs.append(af)

        # 合并文件
        if len(dfs) == 1:
            merged = dfs[0]
        else:
            merged = dfs[0]
            for df in dfs[1:]:
                merged = pd.merge(merged, df, on=merge_cols, how='outer')

        # 排序
        merged = merged.sort_values(by=["ID1", "ID2"], ignore_index=True)

        # 输出文件
        header_note = (
            "MEGADOCK PPI Score > 12, 80% probability of interaction; "
            "> 10, 50% probability of interaction; "
            "> 8, 10% probability of interaction. "
            "HDOCK docking score < -200, high probability of interaction. "
            "pTM + ipTM > 0.75, indicating that the model's predictions of protein-protein interface "
            "interactions and overall complex structure are highly reliable."
        )

        with open(args.output, "w", encoding="utf-8") as f:
            f.write(header_note + "\n")
            merged.to_csv(f, sep="\t", index=False)

        print(f"✅ 合并完成，输出文件：{args.output}")

    except Exception as e:
        print(f"❌ 出错：{e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
