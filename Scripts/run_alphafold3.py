#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import json
import subprocess
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

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

def run_docker_on_json(json_path, output_path, model_dir, db_dir, docker_image="alphafold3", step="Prediction", pdbs_dir=None):
    json_name = os.path.basename(json_path)
    protein_id = json_name.replace("_data.json", "").replace(".json", "")

    print(f"[INFO] Running AlphaFold3 ({step}) for {protein_id}...")

    # 根据步骤添加额外参数
    extra_args = ""
    if step == "Msa":
        extra_args = "--norun_inference"
    elif step == "Inference":
        extra_args = "--norun_data_pipeline"

    # 构建 docker 命令
    cmd = f"""
    docker run --rm {"--gpus all" if step != "Msa" else ""} \\
        -v {os.path.abspath(json_path)}:/root/input.json \\
        -v {os.path.abspath(output_path)}:/root/af_output \\
        -v {os.path.abspath(model_dir)}:/root/models \\
        -v {os.path.abspath(db_dir)}:/root/public_databases \\
        {docker_image} \\
        python run_alphafold.py \\
        --json_path=/root/input.json \\
        --model_dir=/root/models \\
        --output_dir=/root/af_output \\
        {extra_args}
    """
    subprocess.run(cmd, shell=True)

    # 推理步骤才进行 cif → pdb 转换
    if step in ["Inference", "Prediction"]:
        cif_file = os.path.join(output_path, f"{protein_id}/{protein_id}_model.cif")
        if pdbs_dir:
            os.makedirs(pdbs_dir, exist_ok=True)
            pdb_file = os.path.join(pdbs_dir, f"{protein_id}.pdb")

        if os.path.exists(cif_file):
            pymol_cmd = f'pymol -c -q -d "load {cif_file}; save {pdb_file}; quit;"'
            subprocess.run(pymol_cmd, shell=True)
            print(f"[INFO] Converted {protein_id}_model.cif → {protein_id}.pdb")
        else:
            print(f"[WARNING] {cif_file} not found. Skipping PDB conversion.")

def main():
    parser = argparse.ArgumentParser(description="Run AlphaFold3 in MSA/Inference/Prediction mode.")
    parser.add_argument("-s","--step", required=True, choices=["Msa", "Inference", "Prediction"], help="Execution step: Msa, Inference, or Prediction.")
    parser.add_argument("-fa","--fasta", help="Input protein FASTA file (required for Msa and Prediction)")
    parser.add_argument("-j","--json_dir", required=True, help="Directory to store or read JSON files")
    parser.add_argument("-od","--output_dir", required=True, help="Directory to store outputs")
    parser.add_argument("-p","--parameter_dir", required=True, help="AlphaFold3 model parameter directory")
    parser.add_argument("-d","--database_dir", required=True, help="AlphaFold3 public database directory")
    parser.add_argument("-i","--docker_image", default="alphafold3", help="Docker image name")
    parser.add_argument("-n","--num_workers", type=int, default=6, help="Number of concurrent MSA jobs (Msa step only)")

    args = parser.parse_args()

    # 参数验证
    if args.step in ["Msa", "Prediction"] and not args.fasta:
        print("[ERROR] --fasta is required for step Msa or Prediction.")
        sys.exit(1)

    if args.step != "Msa" and any(arg.startswith("--num_workers") for arg in sys.argv):
        print("[ERROR] --num_workers can only be specified with --step Msa.")
        sys.exit(1)

    os.makedirs(args.json_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    pdbs_dir = os.path.join(args.output_dir, "pdbs")
    if args.step in ["Prediction", "Inference"]:
        os.makedirs(pdbs_dir, exist_ok=True)

    # FASTA → JSON（Msa, Prediction）
    if args.step in ["Msa", "Prediction"]:
        sequences = parse_fasta(args.fasta)
        for pid, seq in sequences.items():
            json_obj = convert_to_json_format(pid, seq)
            json_path = os.path.join(args.json_dir, f"{pid}.json")
            with open(json_path, 'w') as f:
                json.dump(json_obj, f, indent=2)

        json_tasks = [
            os.path.join(args.json_dir, f"{pid}.json") for pid in sequences
        ]
    else:
        # Inference 阶段：查找 json_dir/*/*_data.json
        json_tasks = []
        for entry in os.listdir(args.json_dir):
            subdir = os.path.join(args.json_dir, entry)
            if os.path.isdir(subdir):
                data_json_path = os.path.join(subdir, f"{entry}_data.json")
                if os.path.isfile(data_json_path):
                    json_tasks.append(data_json_path)

        if not json_tasks:
            print("[ERROR] No *_data.json files found in output_dir for Inference step.")
            sys.exit(1)

    # 执行任务
    if args.step == "Msa":
        print(f"[INFO] Running MSA with {args.num_workers} concurrent jobs...")
        with ThreadPoolExecutor(max_workers=args.num_workers) as executor:
            futures = []
            for path in json_tasks:
                json_name = os.path.basename(path)
                protein_id = json_name.replace(".json", "")
                data_json_path = os.path.join(args.output_dir, f"{protein_id}/{protein_id}_data.json")

                if os.path.exists(data_json_path):
                    print(f"[SKIP] MSA for {protein_id} already exists. Skipping...")
                    continue

                futures.append(
                    executor.submit(
                        run_docker_on_json,
                        json_path=path,
                        output_path=args.output_dir,
                        model_dir=args.parameter_dir,
                        db_dir=args.database_dir,
                        docker_image=args.docker_image,
                        step=args.step
                    )
                )
            for future in as_completed(futures):
                future.result()
    else:
        for path in json_tasks:
            json_name = os.path.basename(path)
            protein_id = json_name.replace("_data.json", "").replace(".json", "")
            cif_file = os.path.join(args.output_dir, f"{protein_id}/{protein_id}_model.cif")

            if os.path.exists(cif_file):
                print(f"[SKIP] Result for {protein_id} already exists. Skipping...")
                continue

            run_docker_on_json(
                json_path=path,
                output_path=args.output_dir,
                model_dir=args.parameter_dir,
                db_dir=args.database_dir,
                docker_image=args.docker_image,
                step=args.step,
                pdbs_dir=pdbs_dir
            )

    print(f"[DONE] AlphaFold3 step `{args.step}` completed.")

if __name__ == "__main__":
    main()
