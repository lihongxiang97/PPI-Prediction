#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Hongxiang Li, 2025/07/09

# This script runs AlphaFold3 for protein structure prediction.

import os
import argparse
import json
import subprocess

# 1. 解析FASTA文件
def parse_fasta(fasta_file):
    sequences = {}
    with open(fasta_file, "r") as file:
        protein_id = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if protein_id:
                    sequences[protein_id] = "".join(sequence)
                protein_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
        if protein_id:
            sequences[protein_id] = "".join(sequence)
    return sequences

# 2. 构建 AlphaFold3 所需的 JSON 格式
def convert_to_json_format(protein_id, sequence):
    return {
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
        "name": protein_id,
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": sequence
                }
            }
        ]
    }

# 3. 用docker运行AF3 + 用pymol转cif为pdb
def run_docker_on_json(json_path, output_path, model_dir, db_dir, docker_image="alphafold3"):
    json_name = os.path.basename(json_path)
    protein_id = json_name.replace(".json", "")

    print(f"[INFO] Running AlphaFold3 for {protein_id}...")

    cmd = f"""
    docker run --rm \
        --gpus all \
        -v {os.path.abspath(json_path)}:/root/input.json \
        -v {os.path.abspath(output_path)}:/root/af_output \
        -v {os.path.abspath(model_dir)}:/root/models \
        -v {os.path.abspath(db_dir)}:/root/public_databases \
        {docker_image} \
        python run_alphafold.py \
        --json_path=/root/input.json \
        --model_dir=/root/models \
        --output_dir=/root/af_output
    """
    subprocess.run(cmd, shell=True)

    # 自动转换cif为pdb
    cif_file = os.path.join(output_path, f"{protein_id}/{protein_id}_model.cif")
    pdb_file = os.path.join(output_path, f"{protein_id}/{protein_id}.pdb")

    if os.path.exists(cif_file):
        pymol_cmd = f'pymol -c -q -d "load {cif_file}; save {pdb_file}; quit;"'
        subprocess.run(pymol_cmd, shell=True)
        print(f"[INFO] Converted {protein_id}_model.cif → {protein_id}.pdb")
    else:
        print(f"[WARNING] {cif_file} not found. Skipping PDB conversion.")

# 4. 主函数
def main():
    parser = argparse.ArgumentParser(description="Batch run AlphaFold3 from FASTA.")
    parser.add_argument("--fasta", required=True, help="Input protein FASTA file")
    parser.add_argument("--json_dir", required=True, help="Directory to store converted JSON files")
    parser.add_argument("--output_dir", required=True, help="Directory to store predicted PDBs")
    parser.add_argument("--model_dir", required=True, help="AlphaFold3 model parameter directory")
    parser.add_argument("--database_dir", required=True, help="AlphaFold3 public database directory")
    parser.add_argument("--docker_image", default="alphafold3", help="Docker image name")

    args = parser.parse_args()

    os.makedirs(args.json_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    sequences = parse_fasta(args.fasta)
    for pid, seq in sequences.items():
        json_obj = convert_to_json_format(pid, seq)
        json_path = os.path.join(args.json_dir, f"{pid}.json")
        with open(json_path, 'w') as f:
            json.dump(json_obj, f, indent=2)

    for pid in sequences:
        json_file = os.path.join(args.json_dir, f"{pid}.json")
        run_docker_on_json(
            json_path=json_file,
            output_path=args.output_dir,
            model_dir=args.model_dir,
            db_dir=args.database_dir,
            docker_image=args.docker_image
        )

    print("[DONE] All predictions and conversions completed.")

if __name__ == "__main__":
    main()