#!/bin/bash
#$ -l highp,h_data=20G,h_rt=24:00:00
#$ -wd <insert path to working directory>
#$ -o <insert path to log directory>
#$ -e <insert path to log directory>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_CAQU_SNPRelate_VCF2GDS

## Import Packages

source <insert path to miniconda>
conda activate my_SNPRelate

set -o pipefail

## Main

Rscript step04_a_CAQU_SNPRelate_VCF2GDS_20240617.R 

## Clean Up
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate

