#!/bin/bash
#$ -l highp,h_data=20G,h_rt=24:00:00
#$ -wd <insert path to working directory>
#$ -o <insert path to log directory>
#$ -e <insert path to log directory>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_c_CAQU_SNPRelateIBD_KING_20240617

## Import Packages

source <insert path to miniconda>
conda activate my_SNPRelate

set -o pipefail

## Main

Rscript step04_b_CAQU_SNPRelateIBD_KING_20240617.R 

## Clean Up
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate

