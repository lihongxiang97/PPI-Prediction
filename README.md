# PPI-Prediction

A protein-protein interaction (PPI) prediction pipeline using MEGADOCK, HDOCK, and AlphaFold3.

| | |
| --- | --- |
| Authors | Hongxiang Li ([Hongxiang Li](https://github.com/lihongxiang97)) |
| Email   | <qiaoxin@njau.edu.cn> |

## The schematic diagram of DupGen_finder pipeline
<p align="center">
<img src="https://github.com/lihongxiang97/PPI-Prediction/blob/main/data/flowchart.png"  height="500" width="300">
<p align="center">
  
<t align="center">
  Figure 1: The flowchart of PPI-prediction pipeline
</t>

## Contents
* [Dependencies](#dependencies)
* [Installation](#installation)
* [Preparing input files](#preparing-input-files)
* [Running](#running)
* [Result Files](#result-files)
* [Citation](#citation)

## Dependencies

- [MEGADOCK](https://github.com/akiyamalab/MEGADOCK)
- [HDOCK](http://hdock.phys.hust.edu.cn/)
- [AlphaFold3](https://github.com/google-deepmind/alphafold3)

## Installation

```bash
cd ~/software  # or any directory of your choice
git clone git@github.com:lihongxiang97/PPI-Prediction.git
```

**\*\*Note\*\***
Before you use PPI-Prediction, you must install MEGADOCK, HDOCK and Alphafold3.
