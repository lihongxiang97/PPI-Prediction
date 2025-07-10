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

```bash
python scripts/run_alphafold3.py \
    --fasta data/pep.fa \
    --output_dir data/pdbs/
```
This script will generate a .pdb structure for each protein ID in the FASTA file.

The output directory (e.g., data/pdbs/) will contain all predicted structures.

### Step 2: Run MEGADOCK for rigid-body docking
After protein structures are ready, MEGADOCK is used to perform fast rigid-body docking for each protein pair.

```bash
python scripts/run_megadock.py \
    --pair_file data/Protein_pair.list \
    --pdb_dir data/pdbs/ \
    --output_dir results/megadock/
```

### Step 3: Run HDOCK for hybrid docking
HDOCK is then used to provide another set of docking scores using a hybrid algorithm.

```bash
python scripts/run_hdock.py \
    --pair_file data/Protein_pair.list \
    --pdb_dir data/pdbs/ \
    --output_dir results/hdock/
```
💡 Both run_megadock.py and run_hdock.py will automatically prepare the proper input format required by each tool.
### Step 4: Use AlphaFold3 to predict complex structures (optional, slower)
If desired, AlphaFold3 can also be used to directly predict the protein complex structure for each pair. This provides a third type of interaction confidence metric, such as predicted interface pLDDT or pDockQ.

```bash
python scripts/run_alphafold3_complex.py \
    --pair_file data/Protein_pair.list \
    --fasta data/pep.fa \
    --output_dir results/af3_complex/
```
⚠️ This step is computationally expensive and requires GPU.

### Step 5: Merge and filter scores
Finally, the three sources of interaction scores are merged and filtered to generate a list of high-confidence protein interactions.

```bash
python scripts/score_merge.py \
    --megadock results/megadock/scores.tsv \
    --hdock results/hdock/scores.tsv \
    --af3 results/af3_complex/scores.tsv \
    --output results/merged_scores.tsv
```
