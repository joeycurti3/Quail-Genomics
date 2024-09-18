#!/bin/bash
#$ -l highp,h_data=20G,h_rt=6:00:00
#$ -wd <insert working directory path>
#$ -o <insert log directory path>
#$ -e <insert log directory path>
#$ -m bae
#$ -M 1joeynik
#$ -N step04_a_CAQU_SNPRelatePCA_ld2_maf05_unrelated_SAMO

## Import Packages

source <insert path to miniconda>
conda activate my_SNPRelate

set -o pipefail

## Main

Rscript step04_a_CAQU_SNPRelatePCA_ld2_maf05_unrelated_SAMO_20240822.R

## Clean Up
echo -e "[$(date "+%Y-%m-%d %T")] Done"

conda deactivate

