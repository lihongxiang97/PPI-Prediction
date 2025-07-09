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

- docker
- pymol
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

## Step 1: Predict 3D structures of all proteins using AlphaFold3
You need to generate 3D structures of all individual proteins using AlphaFold3. This step must be completed before docking.

```bash
python scripts/run_alphafold3.py \
    --fasta data/pep.fa \
    --output_dir data/pdbs/
```
This script will generate a .pdb structure for each protein ID in the FASTA file.

The output directory (e.g., data/pdbs/) will contain all predicted structures.

## Step 2: Run MEGADOCK for rigid-body docking
After protein structures are ready, MEGADOCK is used to perform fast rigid-body docking for each protein pair.

```bash
python scripts/run_megadock.py \
    --pair_file data/Protein_pair.list \
    --pdb_dir data/pdbs/ \
    --output_dir results/megadock/
```
