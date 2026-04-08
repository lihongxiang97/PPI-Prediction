#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import shutil
import subprocess
import sys


def check_file(path, label):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{label} not found: {path}")


def check_dir(path, label):
    if not os.path.isdir(path):
        raise FileNotFoundError(f"{label} not found: {path}")


def check_command(cmd, label):
    if shutil.which(cmd) is None:
        raise EnvironmentError(f"{label} is not available in PATH: {cmd}")


def maybe_check_hdock(hdock_path):
    if hdock_path:
        hdock_exec = os.path.join(hdock_path, "hdock")
        createpl_exec = os.path.join(hdock_path, "createpl")
        check_file(hdock_exec, "HDOCK executable")
        check_file(createpl_exec, "HDOCK createpl executable")
    else:
        check_command("hdock", "HDOCK executable")
        check_command("createpl", "HDOCK createpl executable")


def build_paths(base_dir):
    return {
        "af_json_dir": os.path.join(base_dir, "af_json"),
        "af_output_dir": os.path.join(base_dir, "af_output"),
        "megadock_output_dir": os.path.join(base_dir, "megadock_out"),
        "megadock_result": os.path.join(base_dir, "megadock.tsv"),
        "hdock_output_dir": os.path.join(base_dir, "hdock_out"),
        "hdock_result": os.path.join(base_dir, "hdock.tsv"),
        "af_complex_json_dir": os.path.join(base_dir, "af_complex_json"),
        "af_complex_output_dir": os.path.join(base_dir, "af_complex_out"),
        "af_complex_result": os.path.join(base_dir, "af_complex.tsv"),
        "merged_result": os.path.join(base_dir, "merged_scores.tsv"),
    }


def run_command(cmd, description):
    print(f"[STEP] {description}")
    print("[CMD] " + " ".join(cmd))
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(f"{description} failed with exit code {result.returncode}")


def validate_environment(args, paths):
    check_file(args.pair_list, "Protein pair list")
    check_file(args.fasta, "Protein FASTA")

    if not args.skip_single_af or not args.skip_complex_af:
        check_dir(args.parameter_dir, "AlphaFold3 parameter directory")
        check_dir(args.database_dir, "AlphaFold3 database directory")
        check_command("docker", "Docker")

    if not args.skip_hdock:
        maybe_check_hdock(args.hdock_path)

    if args.convert_complex_pdb or not args.skip_single_af:
        check_command("pymol", "PyMOL")


def run_single_alphafold(script_path, args, paths):
    cmd = [
        sys.executable,
        script_path,
        "--step", args.af_step,
        "--fasta", args.fasta,
        "--json_dir", paths["af_json_dir"],
        "--output_dir", paths["af_output_dir"],
        "--parameter_dir", args.parameter_dir,
        "--database_dir", args.database_dir,
        "--docker_image", args.af_docker_image,
    ]
    if args.af_step == "Msa":
        cmd.extend(["--num_workers", str(args.num_workers)])
    run_command(cmd, "AlphaFold3 single-protein prediction")


def run_megadock(script_path, args, paths):
    pdb_dir = os.path.join(paths["af_output_dir"], "pdbs")
    check_dir(pdb_dir, "AlphaFold3 single-protein PDB directory")
    cmd = [
        sys.executable,
        script_path,
        "--pair_list", args.pair_list,
        "--pdb_dir", pdb_dir,
        "--output_dir", paths["megadock_output_dir"],
        "--result_file", paths["megadock_result"],
        "--docker_image", args.megadock_docker_image,
        "-N", str(args.megadock_decoys),
        "-t", str(args.megadock_fft_threads),
        "-e", str(args.megadock_cpu_cores),
    ]
    run_command(cmd, "MEGADOCK docking")


def run_hdock(script_path, args, paths):
    pdb_dir = os.path.join(paths["af_output_dir"], "pdbs")
    check_dir(pdb_dir, "AlphaFold3 single-protein PDB directory")
    cmd = [
        sys.executable,
        script_path,
        "--pair_list", args.pair_list,
        "--pdb_dir", pdb_dir,
        "--output_dir", paths["hdock_output_dir"],
        "--result_file", paths["hdock_result"],
        "--threads", str(args.hdock_threads),
    ]
    if args.hdock_path:
        cmd.extend(["--hdock_path", args.hdock_path])
    run_command(cmd, "HDOCK docking")


def run_complex_alphafold(script_path, args, paths):
    cmd = [
        sys.executable,
        script_path,
        "--pair_list", args.pair_list,
        "--fasta", args.fasta,
        "--json_dir", paths["af_complex_json_dir"],
        "--output_dir", paths["af_complex_output_dir"],
        "--model_dir", args.parameter_dir,
        "--database_dir", args.database_dir,
        "--docker_image", args.af_docker_image,
        "--outfile", paths["af_complex_result"],
    ]
    if args.convert_complex_pdb:
        cmd.append("--convert_pdb")
    run_command(cmd, "AlphaFold3 complex prediction")


def run_merge(script_path, args, paths):
    cmd = [
        sys.executable,
        script_path,
        "--output", paths["merged_result"],
    ]
    if not args.skip_megadock and os.path.exists(paths["megadock_result"]):
        cmd.extend(["--megadock", paths["megadock_result"]])
    if not args.skip_hdock and os.path.exists(paths["hdock_result"]):
        cmd.extend(["--hdock", paths["hdock_result"]])
    if not args.skip_complex_af and os.path.exists(paths["af_complex_result"]):
        cmd.extend(["--af", paths["af_complex_result"]])
    run_command(cmd, "Merge score tables")


def parse_args():
    parser = argparse.ArgumentParser(
        description="One-click launcher for the PPI-Prediction pipeline."
    )
    parser.add_argument("-l", "--pair_list", required=True, help="Protein pair list file")
    parser.add_argument("-fa", "--fasta", required=True, help="Protein FASTA file")
    parser.add_argument("-o", "--work_dir", default="pipeline_run", help="Working directory for all outputs")
    parser.add_argument("-p", "--parameter_dir", help="AlphaFold3 model parameter directory")
    parser.add_argument("-d", "--database_dir", help="AlphaFold3 database directory")
    parser.add_argument("--af_step", choices=["Prediction", "Msa"], default="Prediction",
                        help="Single-protein AlphaFold3 mode. Use Prediction for one-step run, Msa for MSA-only.")
    parser.add_argument("--af_docker_image", default="alphafold3", help="AlphaFold3 Docker image name")
    parser.add_argument("--megadock_docker_image", default="hub.rat.dev/akiyamalab/megadock:gpu",
                        help="MEGADOCK Docker image name")
    parser.add_argument("--hdock_path", default=None, help="Optional directory containing hdock and createpl")
    parser.add_argument("--num_workers", type=int, default=6, help="Concurrent MSA jobs for AlphaFold3 Msa mode")
    parser.add_argument("--megadock_decoys", type=int, default=10800, help="MEGADOCK decoy number")
    parser.add_argument("--megadock_fft_threads", type=int, default=3, help="MEGADOCK FFT thread number")
    parser.add_argument("--megadock_cpu_cores", type=int, default=32, help="MEGADOCK OMP_NUM_THREADS value")
    parser.add_argument("--hdock_threads", type=int, default=8, help="Parallel HDOCK jobs")
    parser.add_argument("--convert_complex_pdb", action="store_true",
                        help="Convert AlphaFold3 complex CIF files to PDB with PyMOL")
    parser.add_argument("--skip_single_af", action="store_true",
                        help="Skip single-protein AlphaFold3 and reuse existing work_dir/af_output/pdbs")
    parser.add_argument("--skip_megadock", action="store_true", help="Skip MEGADOCK step")
    parser.add_argument("--skip_hdock", action="store_true", help="Skip HDOCK step")
    parser.add_argument("--skip_complex_af", action="store_true",
                        help="Skip AlphaFold3 complex prediction step")
    parser.add_argument("--skip_merge", action="store_true", help="Skip score merge step")
    return parser.parse_args()


def main():
    args = parse_args()

    args.pair_list = os.path.abspath(args.pair_list)
    args.fasta = os.path.abspath(args.fasta)
    args.work_dir = os.path.abspath(args.work_dir)
    args.parameter_dir = os.path.abspath(args.parameter_dir) if args.parameter_dir else None
    args.database_dir = os.path.abspath(args.database_dir) if args.database_dir else None
    args.hdock_path = os.path.abspath(args.hdock_path) if args.hdock_path else None

    if (not args.skip_single_af or not args.skip_complex_af) and (not args.parameter_dir or not args.database_dir):
        raise ValueError("--parameter_dir and --database_dir are required unless all AlphaFold3 steps are skipped")

    paths = build_paths(args.work_dir)
    os.makedirs(args.work_dir, exist_ok=True)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    scripts = {
        "single_af": os.path.join(script_dir, "run_alphafold3.py"),
        "megadock": os.path.join(script_dir, "run_megadock.py"),
        "hdock": os.path.join(script_dir, "run_hdock.py"),
        "complex_af": os.path.join(script_dir, "run_alphafold3_complex.py"),
        "merge": os.path.join(script_dir, "merge_score.py"),
    }

    try:
        validate_environment(args, paths)

        if not args.skip_single_af:
            run_single_alphafold(scripts["single_af"], args, paths)

        if not args.skip_megadock:
            run_megadock(scripts["megadock"], args, paths)

        if not args.skip_hdock:
            run_hdock(scripts["hdock"], args, paths)

        if not args.skip_complex_af:
            run_complex_alphafold(scripts["complex_af"], args, paths)

        if not args.skip_merge:
            run_merge(scripts["merge"], args, paths)

    except Exception as exc:
        print(f"[ERROR] {exc}")
        sys.exit(1)

    print("[DONE] Pipeline completed successfully.")
    print(f"[INFO] Working directory: {args.work_dir}")
    if not args.skip_merge:
        print(f"[INFO] Merged result: {paths['merged_result']}")


if __name__ == "__main__":
    main()
