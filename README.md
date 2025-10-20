# PPI-Prediction

A protein-protein interaction (PPI) prediction pipeline using MEGADOCK, HDOCK, and AlphaFold3.

| | |
| --- | --- |
| Authors | Hongxiang Li ([Hongxiang Li](https://github.com/lihongxiang97)) |
| Email   | <2022204047@njau.stu.edu.cn> |

## The schematic diagram of DupGen_finder pipeline
<p align="center">
<img src="https://github.com/lihongxiang97/PPI-Prediction/blob/main/data/flowchart.png"  height="800">
<p align="center">
  
<t align="center">
  The flowchart of PPI-prediction pipeline
</t>

## Contents
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Preparing input files](#preparing-input-files)
* [Running](#running)
* [Result Files](#result-files)
* [Citation](#citation)

## Dependencies

- python==3.12.9
- pymol-open-source==3.1.0
- docker
- [MEGADOCK](https://github.com/akiyamalab/MEGADOCK)
- [HDOCK](http://hdock.phys.hust.edu.cn/)
- [AlphaFold3](https://github.com/google-deepmind/alphafold3)

## Installation

```bash
cd ~/software  # or any directory of your choice
git clone https://github.com/lihongxiang97/PPI-Prediction.git
```

**\*\*Note\*\***
Before you use PPI-Prediction, you must install MEGADOCK, HDOCK and Alphafold3.

## Preparing input files

Protein ID pair list and protein sequence fasta file are required for running PPI-prediction successfully.

1.```Protein_pair.list``` is in the following format (tab separated):
```
ID1  ID2
ID1  ID3
ID2  ID4
ID4  ID5
```

2. ```pep.fa``` is in the following format:
```
>ID1
sequence.....
>ID2
sequence.....
>ID3
sequence.....
>ID4
sequence.....
>ID5
sequence.....
```

## Running
The PPI-Prediction pipeline proceeds in the following main steps:

### Step 1: Predict 3D structures of all proteins using AlphaFold3
You need to generate 3D structures of all individual proteins using AlphaFold3. This step must be completed before docking.
The running script is designed based on the Docker installations of AlphaFold3 and MEGADOCK.

```bash
python Scripts/run_alphafold3.py -h
usage: run_alphafold3.py [-h] -s {Msa,Inference,Prediction} [-fa FASTA] -j JSON_DIR -od OUTPUT_DIR -p PARAMETER_DIR -d DATABASE_DIR [-i DOCKER_IMAGE] [-n NUM_WORKERS]

Run AlphaFold3 in MSA/Inference/Prediction mode.

optional arguments:
  -h, --help            show this help message and exit
  -s {Msa,Inference,Prediction}, --step {Msa,Inference,Prediction}
                        Execution step: Msa, Inference, or Prediction.
  -fa FASTA, --fasta FASTA
                        Input protein FASTA file (required for Msa and Prediction)
  -j JSON_DIR, --json_dir JSON_DIR
                        Directory to store or read JSON files
  -od OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to store outputs
  -p PARAMETER_DIR, --parameter_dir PARAMETER_DIR
                        AlphaFold3 model parameter directory
  -d DATABASE_DIR, --database_dir DATABASE_DIR
                        AlphaFold3 public database directory
  -i DOCKER_IMAGE, --docker_image DOCKER_IMAGE
                        Docker image name
  -n NUM_WORKERS, --num_workers NUM_WORKERS
                        Number of concurrent MSA jobs (Msa step only)
```
Sample:
```bash
python scripts/run_alphafold3.py \
    --step Prediction \
    --fasta data/pep.fa \
    --json_dir ./json_dir \
    --output_dir ./af_output \
    --model_dir your/parameters/path \
    --database_dir your/database/path \
    --docker_image your_docker_image_name
```
This script will generate a .pdb structure for each protein ID in the FASTA file.

The output directory (e.g., ./af_output) will contain all predicted structures.

### Step 2: Run MEGADOCK for rigid-body docking
After protein structures are ready, MEGADOCK is used to perform fast rigid-body docking for each protein pair.

```bash
usage: run_megadock.py [-h] -l PAIR_LIST -d PDB_DIR -od OUTPUT_DIR [-r RESULT_FILE] [-i DOCKER_IMAGE] [-N N] [-t T] [-e E]

Run MEGADOCK for PPI prediction

optional arguments:
  -h, --help            show this help message and exit
  -l PAIR_LIST, --pair_list PAIR_LIST
                        Protein pair list file (ID1 ID2)
  -d PDB_DIR, --pdb_dir PDB_DIR
                        Directory containing PDB files
  -od OUTPUT_DIR, --output_dir OUTPUT_DIR
                        MEGADOCK output directory
  -r RESULT_FILE, --result_file RESULT_FILE
                        File to save MEGADOCK scores
  -i DOCKER_IMAGE, --docker_image DOCKER_IMAGE
                        Docker image for MEGADOCK
  -N N                  Number of decoys, default 10800
  -t T                  Thread number for FFT, default 3
  -e E                  Number of CPU cores (OMP_NUM_THREADS), default 32
```
Sample:
```bash
python scripts/run_megadock.py \
    -l data/Protein_pair.list \
    -d af_output/pdbs \
    -od megadock_out \
    -r megadock.tsv \
    -i hub.rat.dev/akiyamalab/megadock:gpu
```

### Step 3: Run HDOCK for hybrid docking
HDOCK is then used to provide another set of docking scores using a hybrid algorithm.

```bash
usage: run_hdock.py [-h] -l PAIR_LIST -d PDB_DIR -od OUTPUT_DIR [-p HDOCK_PATH] [-r RESULT_FILE] [-t THREADS]

Run HDOCK for protein pairs in parallel

optional arguments:
  -h, --help            show this help message and exit
  -l PAIR_LIST, --pair_list PAIR_LIST
                        Protein pair list file
  -d PDB_DIR, --pdb_dir PDB_DIR
                        Directory with PDB files
  -od OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to store HDOCK outputs
  -p HDOCK_PATH, --hdock_path HDOCK_PATH
                        Optional path to HDOCK executables
  -r RESULT_FILE, --result_file RESULT_FILE
                        Output result file
  -t THREADS, --threads THREADS
                        Number of parallel HDOCK tasks (default: 8)
```
Sample:
```bash
python scripts/run_hdock.py \
    --pair_list data/Protein_pair.list \
    --pdb_dir af_output/pdbs \
    --output_dir hdock_out \
    --result_file hdock.tsv \
    -t 10
```
💡 Both run_megadock.py and run_hdock.py will automatically prepare the proper input format required by each tool.
### Step 4: Use AlphaFold3 to predict complex structures (optional, slower)
If desired, AlphaFold3 can also be used to directly predict the protein complex structure for each pair. This provides a third type of interaction confidence metric, such as predicted interface pLDDT or pDockQ.

```bash
usage: run_alphafold3_complex.py [-h] -l PAIR_LIST -fa FASTA -jd JSON_DIR -od OUTPUT_DIR -p MODEL_DIR -d DATABASE_DIR [-i DOCKER_IMAGE] [--convert_pdb] -o OUTFILE

AlphaFold3 Complex Prediction (pair-based)

optional arguments:
  -h, --help            show this help message and exit
  -l PAIR_LIST, --pair_list PAIR_LIST
                        Protein pair list file
  -fa FASTA, --fasta FASTA
                        FASTA file containing all protein sequences
  -jd JSON_DIR, --json_dir JSON_DIR
                        Directory to save generated JSON files
  -od OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory to save prediction results
  -p MODEL_DIR, --model_dir MODEL_DIR
                        AlphaFold3 model parameter directory
  -d DATABASE_DIR, --database_dir DATABASE_DIR
                        AlphaFold3 public database directory
  -i DOCKER_IMAGE, --docker_image DOCKER_IMAGE
                        Docker image name
  --convert_pdb         Convert CIF to PDB using PyMOL
  -o OUTFILE, --outfile OUTFILE
                        Output file to save ptm and iptm results
```
Sample:
```bash
python scripts/run_alphafold3_complex.py \
    -l data/Protein_pair.list \
    -fa data/pep.fa \
    -jd af_c_json \
    -od af_c_out \
    -p your/parameters/path \
    -d your/Database/path \
    --convert_pdb \
    -o af_c.tsv
```
⚠️ This step is computationally expensive and requires GPU.

### Step 5: Merge and filter scores
Finally, the three sources of interaction scores are merged and filtered to generate a list of high-confidence protein interactions.

```bash
python scripts/merge_score.py \
    --megadock megadock.tsv \
    --hdock hdock.tsv \
    --af af_c.tsv \
    --output merged_scores.tsv
```
