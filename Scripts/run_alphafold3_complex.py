#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
import subprocess

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

def convert_complex_to_json(protein_id1, protein_id2, seq1, seq2):
    return {
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1,
        "name": f"{protein_id1}-{protein_id2}",
        "sequences": [
            {
                "protein": {
                    "id": ["A"],
                    "sequence": seq1
                }
            },
            {
                "protein": {
                    "id": ["B"],
                    "sequence": seq2
                }
            }
        ]
    }

def run_docker_prediction(json_path, output_dir, model_dir, db_dir, docker_image, pdbs_dir=None):
    pair_name = os.path.basename(json_path).replace(".json", "")
    lower_pair_name = pair_name.lower()
    pair_out_dir = os.path.join(output_dir, lower_pair_name)
    os.makedirs(pair_out_dir, exist_ok=True)

    cif_file = os.path.join(pair_out_dir, f"{lower_pair_name}_model.cif")
    if os.path.exists(cif_file):
        print(f"[SKIP] Prediction already exists for {pair_name}. Skipping...")
        return

    print(f"[INFO] Predicting complex: {pair_name}...")

    cmd = f"""
    docker run --rm --gpus all \\
        -v {os.path.abspath(json_path)}:/root/input.json \\
        -v {os.path.abspath(output_dir)}:/root/af_output \\
        -v {os.path.abspath(model_dir)}:/root/models \\
        -v {os.path.abspath(db_dir)}:/root/public_databases \\
        {docker_image} \\
        python run_alphafold.py \\
        --json_path=/root/input.json \\
        --model_dir=/root/models \\
        --output_dir=/root/af_output
    """
    subprocess.run(cmd, shell=True)

    if pdbs_dir:
        os.makedirs(pdbs_dir, exist_ok=True)
        pdb_file = os.path.join(pdbs_dir, f"{pair_name}.pdb")
        if os.path.exists(cif_file):
            pymol_cmd = f'pymol -c -q -d "load {cif_file}; save {pdb_file}; quit;"'
            subprocess.run(pymol_cmd, shell=True)
            print(f"[INFO] Converted {pair_name}_model.cif → {pair_name}.pdb")
        else:
            print(f"[WARNING] Missing CIF file for {pair_name}. Cannot convert to PDB.")

def main():
    parser = argparse.ArgumentParser(description="AlphaFold3 Complex Prediction (pair-based)")
    parser.add_argument("-l","--pair_list", required=True, help="Protein pair list file")
    parser.add_argument("-fa","--fasta", required=True, help="FASTA file containing all protein sequences")
    parser.add_argument("-jd","--json_dir", required=True, help="Directory to save generated JSON files")
    parser.add_argument("-od","--output_dir", required=True, help="Directory to save prediction results")
    parser.add_argument("-p","--model_dir", required=True, help="AlphaFold3 model parameter directory")
    parser.add_argument("-d","--database_dir", required=True, help="AlphaFold3 public database directory")
    parser.add_argument("-i","--docker_image", default="alphafold3", help="Docker image name")
    parser.add_argument("--convert_pdb", action="store_true", help="Convert CIF to PDB using PyMOL")

    args = parser.parse_args()

    os.makedirs(args.json_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    pdbs_dir = os.path.join(args.output_dir, "pdbs") if args.convert_pdb else None

    sequences = parse_fasta(args.fasta)

    with open(args.pair_list) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            p1, p2 = parts[0], parts[1]

            if p1 not in sequences or p2 not in sequences:
                print(f"[WARNING] Sequence missing for {p1} or {p2}, skipping...")
                continue

            json_obj = convert_complex_to_json(p1, p2, sequences[p1], sequences[p2])
            json_path = os.path.join(args.json_dir, f"{p1}-{p2}.json")
            with open(json_path, 'w') as fjson:
                json.dump(json_obj, fjson, indent=2)

            run_docker_prediction(
                json_path=json_path,
                output_dir=args.output_dir,
                model_dir=args.model_dir,
                db_dir=args.database_dir,
                docker_image=args.docker_image,
                pdbs_dir=pdbs_dir
            )

    print("[DONE] Complex structure prediction completed.")

if __name__ == "__main__":
    main()
